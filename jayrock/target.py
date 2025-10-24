import warnings

from astropy.time import Time
from astroquery.jplhorizons import Horizons
from jwst_gtvt.jwst_tvt import Ephemeris
import numpy as np
import requests
from rich import print
import rocks

import astro
import jayrock


# Subclass the Ephemeris class to patch the get_moving_target_positions method
# If the package is open for PRs, these changes may be upstreamed in the future
class PatchedEphemeris(Ephemeris):
    def get_moving_target_positions(self, name):
        """Get ephemeris from JPL/HORIZONS.

        Parameters
        ----------
        name : str
            Name of target.

        Returns
        -------
        pd.DataFrame
            DataFrame containing ephemerides at requested dates.

        Notes
        -----
        Changes include
          - id_type to asteroid_name in Horizons query
          - add Vmag, RA, Dec, rh, delta, phase, positional uncertainties to output dataframe
        """
        self.fixed = False

        start = self.start_date.to_value("iso", subfmt="date")
        stop = self.end_date.to_value("iso", subfmt="date")

        obj = Horizons(
            id=name,
            location="500@-170",
            epochs={"start": start, "stop": stop, "step": "1d"},
            id_type="asteroid_name",
        )

        # https://ssd.jpl.nasa.gov/horizons/manual.html#output
        eph = obj.ephemerides(cache=False, quantities="1,9,19,20,24,36")

        self.dataframe["ra"] = eph["RA"].data.data
        self.dataframe["dec"] = eph["DEC"].data.data

        self.dataframe = self.build_dataframe()
        self.dataframe["V"] = eph["V"].data.data
        self.dataframe["rh"] = eph["r"].data.data
        self.dataframe["delta"] = eph["delta"].data.data
        self.dataframe["phase"] = eph["alpha"].data.data
        self.dataframe["ra_3sigma"] = eph["RA_3sigma"].data.data
        self.dataframe["dec_3sigma"] = eph["DEC_3sigma"].data.data

        return self.dataframe


class Target:
    """Define target of observation."""

    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
            Name or desingation of the target to observe. Passed to rocks.id for identification.
        """
        self.rock = rocks.Rock(name)

        if not self.rock.is_valid:
            jayrock.logging.logger.warning(
                "Functionality is limited if target is unknown."
            )
            self.name = name
        else:
            self.name = self.rock.name

    def __repr__(self):
        return f"Target(name={self.name})"

    def compute_ephemeris(self, cycle=None, date_start=None, date_end=None):
        """Compute visibility windows for the target using jwst mtvt and JPL Horizons.

        Parameters
        ----------
        cycle : int, optional
            JWST observation cycle to determine start and end dates of observation.
            Must be provided if date_start and date_end are not given. Default is None.
        date_start : str, optional
            Start date of the observation in ISO format (YYYY-MM-DD). Can be
            omitted if cycle is given. Default is None.
        date_end : str, optional
            End date of the observation in ISO format (YYYY-MM-DD). Can be
            omitted if cycle is given. Default is None.

        Notes
        -----
        Populates the following attributes:
            - visibility : pd.DataFrame
                Visibility windows for the target.
        """
        if cycle is None and date_start is None and date_end is None:
            raise ValueError("Either cycle or date_start and date_end must be given.")

        # Set dates
        if cycle is not None:
            date_start, date_end = jayrock.get_cycle_dates(cycle)

        date_start = Time(date_start, format="iso", scale="utc")
        date_end = Time(date_end, format="iso", scale="utc")

        # suppress ugly warning triggered by mtvt
        # instantiating astropy.time.Time for year 2030
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            eph = PatchedEphemeris(start_date=date_start, end_date=date_end)

        eph.get_moving_target_positions(self.name)

        # Only use dates when target is in field of regard
        vis = eph.dataframe.loc[eph.dataframe["in_FOR"]].copy().reset_index(drop=True)

        # ------
        # Add thermal flux
        target_neatm = self.neatm()

        # Thermal flux at 15 micron in mJy
        vis["thermal_flux"] = vis.apply(
            lambda row: target_neatm.fluxd(
                wave=15,
                geometry=dict(rh=row["rh"], delta=row["delta"], phase=row["phase"]),
            ).value,
            axis=1,
        )

        # ------
        # Add convenience columns and attributes

        # Add observation window column
        vis["window"] = (vis["MJD"].diff() > 1).cumsum() + 1

        # Add date_obs column in iso format
        vis["date_obs"] = vis["MJD"].apply(
            lambda x: Time(x, format="mjd", scale="utc").iso.split()[0]
        )

        self.visibility = vis

        self.vmag_min = self.visibility.V.min()
        self.vmag_max = self.visibility.V.max()

        self.is_visible = len(self.visibility) > 0

        self.n_days_visible = len(self.visibility)
        self.n_windows = self.visibility["window"].max() if self.is_visible else 0

        self.dates_vis_start = [
            window["date_obs"].min() for _, window in self.visibility.groupby("window")
        ]
        self.dates_vis_end = [
            window["date_obs"].max() for _, window in self.visibility.groupby("window")
        ]

    def is_visible_on(self, date_obs):
        """Check if target is visible on a given date.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD).

        Returns
        -------
        bool
            True if target is visible on the given date, False otherwise.
        """
        if not hasattr(self, "visibility"):
            raise AttributeError("Run query_visibility() before is_visible_on().")

        if not any(self.visibility["date_obs"].str.startswith(date_obs)):
            return False
        return True

    def neatm(self, diameter=None, albedo=None, beaming=1.0, emissivity=0.9):
        """Build NEATM model for the target.

        Parameters
        ----------
        diameter : float, optional
            Diameter of the target in km. If None, uses value from rocks database or assumes 30 km.
            Default is None.
        albedo : float, optional
            Albedo of the target. If None, uses value from rocks database or assumes 0.1. Default is None.
        beaming : float, optional
            Beaming parameter. Default is 1.0.
        emissivity : float, optional
            Emissivity of the target. Default is 0.9.

        Returns
        -------
        NEATM
            NEATM model for the target.
        """

        if albedo is None:
            albedo = self.rock.albedo.value if self.rock.albedo else 0.1
            if not self.rock.albedo:
                jayrock.logging.logger.warning(
                    f"Albedo unknown for {self.name}, assuming {albedo}."
                )

        if diameter is None:
            diameter = self.rock.diameter.value if self.rock.diameter else 30.0  # km
            if not self.rock.diameter:
                jayrock.logging.logger.warning(
                    f"Diameter unknown for {self.name}, assuming {diameter} km."
                )

        target_neatm = jayrock.neatm.NEATM(
            diameter=diameter,
            albedo=albedo,
            beaming=1.0,
            emissivity=0.9,
        )
        return target_neatm

    def print_visibility(self):
        """Print visibility windows for the target."""
        print(
            f"\n{self.name}:  Visible for {self.n_days_visible} days during {self.n_windows} window{'' if self.n_windows < 2 else 's'}"
        )

        if "thermal_flux" not in self.visibility.columns and hasattr(
            self, "thermal_flux"
        ):
            self.visibility["thermal_flux"] = self.thermal_flux

        for n, window in self.visibility.groupby("window"):
            window = window.sort_values("date_obs")

            print(
                f"\n  Window {n}: {window.date_obs.min()} -> {window.date_obs.max()} [{len(window)} days]"
            )
            print(
                f"      Vmag:    {window['V'].values[0]:.2f} -> {window['V'].values[-1]:.2f}"
            )
            print(
                f"  Thermal: {window['thermal_flux'].values[0]:.2f} -> {window['thermal_flux'].values[-1]:.2f} mJy [15micron]"
            )

        print(
            f"\nPositional uncertainty: {np.nanmean(self.visibility['ra_3sigma'].values):.3f} arcsec in RA / {np.nanmean(self.visibility['dec_3sigma'].values):.3f} arcsec in DEC\n"
        )

    def _build_config_source_reflected(self, date_obs):
        """Configure the reflected source of the target.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD).

        Returns
        -------
        dict
            Source configuration dictionary to be passed to the pandeia.
        """
        idx_date_obs = self.visibility.date_obs.tolist().index(date_obs)
        self.vmag_date_obs = self.visibility.V.values[idx_date_obs]

        # TODO: Make this configurable
        reflected = {
            "shape": {"geometry": "point"},
            "position": {
                "orientation": 0.0,
                "y_offset": 0.0,
                "x_offset": 0.0,
            },
            "spectrum": {
                "extinction": {
                    "law": "mw_rv_31",
                    "bandpass": "j",
                    "value": 0.0,
                    "unit": "mag",
                },
                "sed": {
                    "key": "g2v",
                    "sed_type": "phoenix",
                    "teff": 5800,
                    "log_g": 4.5,
                    "metallicity": 0.0,
                },
                "extinction_first": True,
                "redshift": 0.0,
                "normalization": {
                    "bandpass": "johnson,v",
                    "norm_flux": self.vmag_date_obs,
                    "norm_fluxunit": "vegamag",
                    "type": "photsys",
                },
                "lines": [],
                "name": "generic source",
            },
            "id": 1,
        }
        return reflected

    def compute_reflectance_spectrum(self, date_obs=None):
        """Compute the reflectance spectrum source at a given date of observation.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD). If None, uses the
            date of brightest apparent V magnitude. Default is None.

        Returns
        -------
        spectrum : AstroSpectrum
            Reflectance spectrum of the target at the given date.
        """
        from pandeia.engine.astro_spectrum import AstroSpectrum
        from pandeia.engine.source import Source

        # TODO: Make this a function
        if date_obs is None:
            if not hasattr(self, "visibility"):
                raise AttributeError(
                    "Run query_visibility() before build_reflectance_spectrum()."
                )
            idx_brightest = self.visibility.V.idxmin()
            date_obs = self.visibility.date_obs.values[idx_brightest]

        config_reflected = self._build_config_source_reflected(date_obs)
        source_reflected = Source(telescope="jwst", config=config_reflected)
        return AstroSpectrum(source_reflected)

    def _build_config_source_thermal(self, date_obs):
        """Configure the thermal source of the target.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD).

        Returns
        -------
        dict
            Source configuration dictionary to be passed to the pandeia.
        """
        # idx_date_obs = self.visibility.date_obs.tolist().index(date_obs)
        # self.thermal_flux_date_obs = self.thermal_flux[idx_date_obs]

        wave, flux = self.compute_thermal_spectrum(date_obs)

        thermal = {
            "shape": {"geometry": "point"},
            "position": {
                "orientation": 0.0,
                "y_offset": 0.0,
                "x_offset": 0.0,
            },
            "spectrum": {
                "extinction": {
                    "law": "mw_rv_31",
                    "bandpass": "j",
                    "value": 0.0,
                    "unit": "mag",
                },
                "sed": {
                    # "sed_type": "blackbody",
                    # "temp": int(self.surface_temperature),
                    "sed_type": "input",
                    "spectrum": [wave, flux],
                },
                "extinction_first": True,
                "redshift": 0.0,
                "normalization": {
                    "type": "none",
                    # "norm_flux": float(self.thermal_flux_date_obs),
                    # "norm_fluxunit": "mjy",
                    # "norm_wave": 15,
                    # "norm_waveunit": "microns",
                    # "type": "at_lambda",
                },
                "lines": [],
                "name": "generic source",
            },
            "id": 2,
        }
        return thermal

    def compute_thermal_spectrum(self, date_obs=None):
        """Build the thermal spectrum source at a given date of observation.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD). If None, uses the
            date of highest thermal flux. Default is None.

        Returns
        -------
        spectrum : AstroSpectrum
            Thermal spectrum of the target at the given date.
        """
        if date_obs is None:
            if not hasattr(self, "visibility"):
                raise AttributeError(
                    "Run query_visibility() before build_thermal_spectrum()."
                )
            idx_highest = self.visibility.thermal_flux_15um_mjy.idxmax()
            date_obs = self.visibility.date_obs.values[idx_highest]
        idx_date_obs = self.visibility.date_obs.tolist().index(date_obs)
        geom = dict(
            rh=self.visibility.rh.values[idx_date_obs],
            delta=self.visibility.delta.values[idx_date_obs],
            phase=self.visibility.phase.values[idx_date_obs],
        )
        target_neatm = self.neatm()
        wave = np.linspace(5.0, 28.0, 200)
        flux = target_neatm.fluxd(wave, geom)
        return wave, flux

    def build_scene(self, date_obs):
        """Build scene for ETC.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD).

        Returns
        -------
        scene : dict
            Scene configuration dictionary to be passed to the pandeia ETC calculation.
        """

        if not self.is_visible_on(date_obs):
            raise ValueError(f"{self.target} is not visible on {date_obs}.")

        if not hasattr(self, "visibility"):
            raise AttributeError("Run query_visibility() before build_scene().")

        if not hasattr(self, "thermal_flux"):
            raise AttributeError("Run query_thermal_flux() before build_scene().")

        reflected = _build_config_source_reflected(date_obs)
        thermal = _build_config_source_thermal(date_obs)
        return [reflected, thermal]

    def plot_spectrum(self, reflected=True, thermal=True, date_obs=None):
        """Plot the reflectance and/or thermal spectrum of the target.

        Parameters
        ----------
        reflected : bool, optional
            Whether to plot the reflectance spectrum. Default is True.
        thermal : bool, optional
            Whether to plot the thermal spectrum. Default is True.
        date_obs : str, optional
            Date of observation in ISO format (YYYY-MM-DD). If None, uses the
            date of brightest apparent V magnitude. Default is None.
        """
        import matplotlib.pyplot as plt

        plt.style.use("default")

        if date_obs is None:
            if not hasattr(self, "visibility"):
                raise AttributeError("Run query_visibility() before plot_spectrum().")
            idx_brightest = self.visibility.V.idxmin()
            date_obs = self.visibility.date_obs.values[idx_brightest]
            print(
                f"Using date of brightest Vmag: {date_obs}, Vmag={self.visibility.V.min():.2f}"
            )

        if reflected:
            spec_reflected = self.compute_reflectance_spectrum(date_obs)
            plt.plot(
                spec_reflected.wave,
                spec_reflected.flux,
                label="Reflected",
                color="blue",
            )

        if thermal:
            spec_thermal = self.compute_thermal_spectrum(date_obs)
            plt.plot(
                spec_thermal.wave,
                spec_thermal.flux,
                label="Thermal",
                color="red",
            )

        plt.xlabel("Wavelength (micron)")
        plt.ylabel("Flux (mJy)")
        plt.title(f"{self.name} Spectrum on {date_obs}")
        plt.legend()
        plt.grid()
        plt.show()
