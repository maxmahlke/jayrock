import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline
import pandas as pd
from pandeia.engine.calc_utils import build_default_calc

import jayrock


class Instrument:
    """Define instrument of observation."""

    def __init__(self, instrument, mode, config=None):
        """
        Parameters
        ----------
        instrument : str
            Instrument name (e.g., 'nirspec', 'miri').
        mode : str
            Observation mode (e.g., 'ifu', 'slit').
        config : dict, optional
            Configuration dictionary for the instrument. If None, a default configuration will be used.
            Configuration dictionary should be the 'configuration' level of the pandeia calculation config.
        """
        self.instrument = instrument.lower()
        self.mode = mode.lower()

        self.config = (
            config
            if config is not None
            else build_default_calc(
                telescope="jwst", instrument=self.instrument, mode=self.mode
            )["configuration"]
        )

    def __repr__(self):
        return f"Instrument(instrument={self.instrument}, mode={self.mode})"

    def print_config(self):
        """Print the observation configuration."""
        print(json.dumps(self.config, indent=4))

    def estimate_ngroups_for_snr(self, snr, target):
        """Estimate the number of groups per integration needed to reach a given SNR for a target.

        Parameters
        ----------
        snr : float
            Desired signal-to-noise ratio.
        target : Target
            Target to observe

        Notes
        -----
        Estimation is based on pre-computed SNR for a sample asteroid.
        """
        target_snr = snr

        PATH_GRID = (
            Path(__file__).parent.parent
            / "data"
            / f"{self.instrument}_{self.mode}_snr_grid.csv"
        )

        if not PATH_GRID.exists():
            raise NotImplementedError(
                f"SNR grid for {self.instrument}/{self.mode} not yet available."
            )

        # The search axes depend on the instrument
        if self.instrument == "nirspec":
            flux, flux_target = "mag", target.vmag_min

        elif self.instrument == "miri":
            flux, flux_target = "logflux", np.log10(target.thermal_flux_max)

        # Load the grid and interpolate
        snr_grid = pd.read_csv(PATH_GRID)
        snr_grid = snr_grid.dropna(subset=["snr_mean"])

        if self.instrument == "miri":
            snr_grid = snr_grid[
                (snr_grid.channel == int(self.config["instrument"]["aperture"][-1]))
                & (snr_grid.disperser == self.config["instrument"]["disperser"])
            ]
        elif self.instrument == "nirspec":
            snr_grid = snr_grid[
                (snr_grid.filter_ == self.config["instrument"]["filter"])
                & (snr_grid.disperser == self.config["instrument"]["disperser"])
            ]

        fluxes = np.sort(snr_grid[flux].unique())
        ngroups = np.sort(snr_grid["ngroup"].unique())

        # Reference calculations were done with nexp = 4
        scale = np.sqrt(self.config["detector"]["nexp"] / 4)
        target_snr /= scale  # TODO: Verify this

        # Pivot the table to create the 2D SNR grid
        snr_grid = snr_grid.pivot(
            index="ngroup", columns=flux, values="snr_mean"
        ).values

        # Interpolate 2d grid to find ngroups for given snr and vmag
        if flux_target < min(fluxes) or flux_target > max(fluxes):
            jayrock.logging.logger.error(
                f"{flux_target} is outside the grid range [{min(fluxes)}, {max(fluxes)}]. Cannot estimate ngroups."
            )

        # Create a 2D spline interpolator for the SNR grid
        interpolator = RectBivariateSpline(ngroups, fluxes, snr_grid, kx=1, ky=1)

        # For a fixed target_flux, find the range of SNRs possible within the grid's exposure times
        min_snr_at_flux = interpolator(ngroups[0], flux_target)[0][0]
        max_snr_at_flux = interpolator(ngroups[-1], flux_target)[0][0]

        if min_snr_at_flux >= target_snr:
            jayrock.logging.logger.info(
                f"Target SNR {target_snr} is met at minimum ngroups = {ngroups[0]}"
            )
            return int(np.ceil(ngroups[0]))

        if target_snr > max_snr_at_flux:
            jayrock.logging.logger.warning(
                f"Target SNR {target_snr} is not achievable for flux {np.power(10,flux_target):.2f}. "
                f"Possible SNR range is [{min_snr_at_flux:.2f}, {max_snr_at_flux:.2f}]. Setting ngroups=10"
            )
            return 10

        snrs_interp = interpolator(ngroups, flux_target).flatten()
        ans = int(np.ceil(np.interp(target_snr, snrs_interp, ngroups)))
        jayrock.logging.logger.info(
            f"Target SNR {target_snr} is met at minimum ngroups = {ngroups[0]}"
        )
        return ans

    @property
    def reference_wavelength(self):
        """Return the reference wavelength of the instrument in microns."""

        if self.instrument == "nirspec":
            return {
                "f070lp": 1.11,
                "f100lp": 1.43,
                "f170lp": 2.415,
                "f290lp": 4.07,
                "CLEAR": 3.5,
            }[self.config["instrument"]["filter"]]

        if self.instrument == "miri":
            return {
                "ch1": {"short": 5.335, "medium": 6.16, "long": 7.109},
                "ch2": {"short": 8.154, "medium": 9.415, "long": 10.85},
                "ch3": {"short": 12.504, "medium": 14.5, "long": 16.745},
                "ch4": {"short": 19.29, "medium": 22.485, "long": 26.22},
            }[self.config["instrument"]["aperture"]][
                self.config["instrument"]["disperser"]
            ]
