import warnings

from astropy.time import Time
from astroquery.jplhorizons import Horizons
from jwst_gtvt.jwst_tvt import Ephemeris
import numpy as np
from pandeia.engine.astro_spectrum import AstroSpectrum
from pandeia.engine.source import Source
from rich import print
import rocks

import jayrock

WAVE_GRID = np.arange(0.45, 30, 0.01)


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

        # If literature values do not exist, populate with defaults
        if not self.rock.albedo:
            self.rock.albedo.value = 0.1
            jayrock.logging.logger.warning(
                f"No albedo value on record for {self.name}. Using default of 0.1."
            )

        if not self.rock.diameter:
            self.rock.diameter.value = 30
            jayrock.logging.logger.warning(
                f"No diameter value on record for {self.name}. Using default of 30."
            )

        self.albedo = self.rock.albedo.value
        self.diameter = self.rock.diameter.value
        self.beaming = 1
        self.emissivity = 0.9

    def __repr__(self):
        return f"Target(name={self.name})"

    def compute_ephemeris(self, cycle=None, date_start=None, date_end=None):
        """Compute target ephemeris using jwst mtvt and JPL Horizons.

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
            - ephemeris : pd.DataFrame
        """
        if cycle is None and date_start is None and date_end is None:
            raise ValueError("Either cycle or date_start and date_end must be given.")

        # Set dates
        if cycle is not None:
            self.date_start, self.date_end = jayrock.get_cycle_dates(cycle)
        else:
            self.date_start = date_start
            self.date_end = date_end

        date_start = Time(self.date_start, format="iso", scale="utc")
        date_end = Time(self.date_end, format="iso", scale="utc")

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
        target_neatm = self.compute_neatm()

        # Thermal flux at 15 micron in mJy
        vis["thermal"] = vis.apply(
            lambda row: target_neatm.fluxd(
                wave=15,
                geometry=dict(rh=row["rh"], delta=row["delta"], phase=row["phase"]),
            ),
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

        self.ephemeris = vis[
            [
                "date_obs",
                "ra",
                "dec",
                "V",
                "rh",
                "delta",
                "phase",
                "ra_3sigma",
                "dec_3sigma",
                "thermal",
                "window",
            ]
        ]

        self.vmag_min = self.ephemeris.V.min()
        self.vmag_max = self.ephemeris.V.max()

        self.is_visible = len(self.ephemeris) > 0

        self.n_days_visible = len(self.ephemeris)
        self.n_windows = self.ephemeris["window"].max() if self.is_visible else 0

        self.dates_vis_start = [
            window["date_obs"].min() for _, window in self.ephemeris.groupby("window")
        ]
        self.dates_vis_end = [
            window["date_obs"].max() for _, window in self.ephemeris.groupby("window")
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
        if not hasattr(self, "ephemeris"):
            raise AttributeError("Run compute_ephemeris() before is_visible_on().")

        if not any(self.ephemeris["date_obs"].str.startswith(date_obs)):
            return False
        return True

    def compute_neatm(self):
        """Build NEATM model for the target.

        Returns
        -------
        NEATM
            NEATM model for the target.
        """
        return jayrock.neatm.NEATM(
            diameter=self.diameter,
            albedo=self.albedo,
            beaming=self.beaming,
            emissivity=self.emissivity,
        )

    def print_ephemeris(self):
        """Print ephemeris of the target."""
        from rich.tree import Tree

        tree = Tree(
            f"({self.rock.number}) {self.name}: Ephemeris from {self.date_start} to {self.date_end}",
        )

        for n, window in self.ephemeris.groupby("window"):
            window = window.sort_values("date_obs")

            branch = tree.add(
                f"[bold]Window {n}:[/bold] {window.date_obs.min()} -> {window.date_obs.max()}"
            )
            branch.add(f"{'Duration':<17}{len(window)} days")
            branch.add(
                f"{'Vmag':<17}{window['V'].values[0]:.2f} -> {window['V'].values[-1]:.2f}"
            )
            branch.add(
                f"{'Thermal @ 15um':<17}{window['thermal'].values[0]:.2f} -> {window['thermal'].values[-1]:.2f} mJy"
            )
        tree.add(
            f"errRA/errDec in arcsec: {np.nanmean(self.ephemeris['ra_3sigma'].values):.3f} / {np.nanmean(self.ephemeris['dec_3sigma'].values):.3f}"
        )
        print(tree)

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
        self.vmag_date_obs = self.ephemeris.loc[
            self.ephemeris.date_obs == date_obs, "V"
        ].values[0]

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

    def compute_reflectance_spectrum(self, date_obs):
        """Compute the reflectance spectrum source at a given date of observation.

        Parameters
        ----------
        date_obs : str
            Date of observation in ISO format (YYYY-MM-DD).

        Returns
        -------
        np.ndarray, np.ndarray
            Wavelengths (micron) and fluxes (mJy) of the reflectance spectrum.
        """
        config_reflected = self._build_config_source_reflected(date_obs)
        source_reflected = Source(telescope="jwst", config=config_reflected)
        spec = AstroSpectrum(source_reflected)
        wave, flux = spec.wave, spec.flux

        # interpolate to WAVE_GRID
        flux = np.interp(WAVE_GRID, wave, flux)
        return WAVE_GRID, flux

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
        self.thermal_date_obs = self.ephemeris.loc[
            self.ephemeris.date_obs == date_obs, "thermal"
        ].values[0]

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
                    # "norm_flux": float(self.thermal),
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
        np.ndarray, np.ndarray
            Wavelengths (micron) and fluxes (mJy) of the thermal spectrum.
        """
        if date_obs is None:
            if not hasattr(self, "ephemeris"):
                raise AttributeError(
                    "Run compute_ephemeris() before build_thermal_spectrum()."
                )
            idx_highest = self.ephemeris.thermal.idxmax()
            date_obs = self.ephemeris.date_obs.values[idx_highest]
        idx_date_obs = self.ephemeris.date_obs.tolist().index(date_obs)
        geom = dict(
            rh=self.ephemeris.rh.values[idx_date_obs],
            delta=self.ephemeris.delta.values[idx_date_obs],
            phase=self.ephemeris.phase.values[idx_date_obs],
        )
        neatm = self.compute_neatm()
        flux = neatm.fluxd(WAVE_GRID, geom)
        return WAVE_GRID, flux

    def get_date_obs(self, at=None):
        """Get list of possible observation dates. Optionally select a date
        based on a condition.

        Parameters
        ----------
        at : str or list of str, optional
           Condition to select date_obs at.
           Must be one of ['vmag_min', 'vmag_max', 'thermal_min', 'thermal_max'].

        Returns
        -------
        list of str or str
            List of possible observation dates. If at is a str or list of
            length 1, returns a single str. If at is a list of length >1,
            returns a list.
        """
        if not hasattr(self, "ephemeris"):
            raise AttributeError("Run compute_ephemeris() before get_date_obs().")

        if at is not None:
            if not isinstance(at, list):
                at = [at]

            date_obs = []

            for at_ in at:
                VALID_AT = ["vmag_min", "vmag_max", "thermal_max", "thermal_min"]

                if at_ not in VALID_AT:
                    raise ValueError(
                        f"Invalid value for at: {at_}. Choose from {VALID_AT}."
                    )

                prop, cond = at_.split("_")
                prop = "V" if prop == "vmag" else "thermal"

                if cond == "min":
                    idx = self.ephemeris[prop].idxmin()
                elif cond == "max":
                    idx = self.ephemeris[prop].idxmax()

                date_obs.append(self.ephemeris.date_obs.values[idx])
            if len(date_obs) == 1:
                return date_obs[0]
            return date_obs
        return self.ephemeris.date_obs.tolist()

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

        if not hasattr(self, "ephemeris"):
            raise AttributeError("Run compute_ephemeris() before build_scene().")

        reflected = self._build_config_source_reflected(date_obs)
        thermal = self._build_config_source_thermal(date_obs)
        return [reflected, thermal]

    def export_spectrum(self, filename, date_obs=None, thermal=True, reflected=True):
        """Export the reflectance and/or thermal spectrum of the target to a file.

        Parameters
        ----------
        filename : str
            Output filename.
        date_obs : str, optional
            Date of observation in ISO format (YYYY-MM-DD). If None, uses the
            date of brightest apparent V magnitude. Default is None.
        thermal : bool, optional
            Whether to export the thermal spectrum. Default is True.
        reflected : bool, optional
            Whether to export the reflectance spectrum. Default is True.
        """

        if not hasattr(self, "ephemeris"):
            raise AttributeError("Run compute_ephemeris() before export_spectrum().")

        wave, flux = WAVE_GRID, np.zeros_like(WAVE_GRID)

        if reflected:
            _, reflected_flux = self.compute_reflectance_spectrum(date_obs=date_obs)
            flux += reflected_flux

        if thermal:
            _, thermal_flux = self.compute_thermal_spectrum(date_obs=date_obs)
            flux += thermal_flux

        # Store to filename
        np.savetxt(
            filename,
            np.column_stack((wave, flux)),
            header="Wavelength(micron) Flux(mJy)",
            delimiter=" ",
            fmt="%.6f",
        )

    def get_vmag(self, date_obs):
        """Get apparent V mag on date of observation."""
        return self.ephemeris.loc[self.ephemeris.date_obs == date_obs, "V"].values[0]

    def get_thermal(self, date_obs, wave=15):
        """Get thermal flux in mJy at wavelength in micron on date of observation."""
        thermal_wave, thermal_flux = self.compute_thermal_spectrum(date_obs=date_obs)
        return np.interp(wave, thermal_wave, thermal_flux)

    def plot_spectrum(self, **kwargs):
        """Plot the spectrum of the target.

        Parameters
        ----------
        target : Target
            Target to plot the spectrum for.
        date_obs : str or list of str
            Date of observation in ISO format (YYYY-MM-DD).
        reflected : bool, optional
            Whether to plot the reflectance spectrum. Default is True.
        thermal : bool, optional
            Whether to plot the thermal spectrum. Default is True.
        """

        if not hasattr(self, "ephemeris"):
            raise AttributeError("Run compute_ephemeris() before plot_spectrum().")

        jayrock.plotting.plot_spectrum(self, **kwargs)
