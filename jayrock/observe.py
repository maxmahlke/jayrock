import collections.abc
import json
from copy import deepcopy
import warnings

import numpy as np
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation

import jayrock


def update(d, u):
    """Recursively update a dictionary. Used to update config dictionaries."""
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


class Observation:
    """JWST observation of a given target with a given instrument."""

    def __init__(self, target, instrument, config, report, date_obs):
        """Configure the observation.

        Parameters
        ----------
        target : Target
            Target to observe.
        instrument : Instrument
            Instrument to use for the observation.
        config : dict
            Configuration dictionary used for the observation.
        report : Report
            Report of the observation.
        date_obs : str
            Date of the observation in ISO format (YYYY-MM-DD).
        """
        self.target = target
        self.instrument = deepcopy(instrument)  # required to protect instrument config
        self.config = config
        self.report = report
        self.date_obs = date_obs

        self.texp = report["scalar"]["total_exposure_time"]
        self.snr = report["scalar"]["sn"]

        self.vmag = target.vmag_date_obs
        self.thermal = target.thermal_date_obs

        self.wave, self.snr_1d = report["1d"]["sn"]
        self.n_partial_saturated = report["1d"]["n_partial_saturated"][1]
        self.n_full_saturated = report["1d"]["n_full_saturated"][1]

        self.partially_saturated = np.array(self.n_partial_saturated) > 0
        self.fully_saturated = np.array(self.n_full_saturated) > 0

        self.is_partially_saturated = np.any(self.n_partial_saturated > 0)
        self.is_fully_saturated = np.any(self.n_full_saturated > 0)
        self.is_saturated = self.is_partially_saturated or self.is_fully_saturated

        if report["warnings"]:
            jayrock.logging.logger.warning(
                f"Observation warnings: {report['warnings']}"
            )

    def print_config(self):
        """Print the observation configuration."""
        print(json.dumps(self.config, indent=4))

    def write_config(self, filename):
        """Write the observation configuration to a JSON file."""
        with open(filename, "w") as f:
            json.dump(config, f, indent=4)

    def plot_snr(self, **kwargs):
        """Plot the SNR as a function of wavelength."""
        return jayrock.plotting.plot_snr(self, **kwargs)


def observe(
    target, instrument, date_obs, scene=None, config_override=None, verbose=True
):
    """Observe a target with a given instrument.

    Parameters
    ----------
    target : Target
        Target to observe.
    instrument : Instrument
        Instrument to use for the observation.
    date_obs : str
        Date of observation in ISO format (YYYY-MM-DD). Must be in the target's
        ephemeris table.
    scene : dict, optional
        Scene dictionary to use for the observation. If None, the scene is built
        from the target properties.
    config_override : dict, optional
        Configuration dictionary to override the default configuration.
    verbose : bool, optional
        Whether to print information about the observation. Default is True.

    Returns
    -------
    Report
        Report of the observation.
    """
    if not hasattr(target, "ephemeris"):
        raise ValueError(
            "Target must have ephemeris information. Run target.compute_ephemeris() first."
        )

    if not target.is_visible_on(date_obs):
        raise ValueError(f"{target} is not visible on {date_obs}. Stopping here.")

    if verbose:
        jayrock.logging.logger.info(
            f"Observing {target} with {instrument.name}|{instrument.mode} on {date_obs}"
        )
        jayrock.logging.logger.info(
            f"{instrument.primary}|{instrument.secondary} - ngroup={instrument.detector.ngroup}|nint={instrument.detector.nint}|nexp={instrument.detector.nexp} - readout={instrument.detector.readout_pattern}"
        )

    # ------
    # Build config
    config = build_default_calc(
        telescope="jwst", instrument=instrument.instrument, mode=instrument.mode
    )
    config["configuration"].update(instrument.config)

    # Build scence from target properties
    config["scene"] = scene if scene is not None else target.build_scene(date_obs)

    # Strategy and background for SNR calculation
    config["strategy"]["reference_wavelength"] = instrument.detector.midpoint
    config["background_level"] = "none"

    if config_override is not None:
        config = update(config, config_override)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = perform_calculation(config, webapp=False)

        snr = report["scalar"]["sn"]
        texp = report["scalar"]["total_exposure_time"]

        if verbose:
            jayrock.logging.logger.info(
                f"SNR={snr:.1f} at {config['strategy']['reference_wavelength']:.2f}μm in {texp / 60:.1f}min"
            )

    return Observation(target, instrument, config, report, date_obs)


def get_cycle_dates(cycle):
    """Get start and end dates of a JWST observation cycle.

    Parameters
    ----------
    cycle : int
        JWST observation cycle number.

    Returns
    -------
    tuple
        Start and end dates of the cycle in ISO format (YYYY-MM-DD).
    """
    CYCLE_DATES = {
        5: ("2026-07-01", "2027-06-30"),
        6: ("2027-07-01", "2028-06-30"),
        7: ("2028-07-01", "2029-06-30"),
        8: ("2029-07-01", "2030-06-30"),
        9: ("2030-07-01", "2031-06-30"),
        10: ("2031-07-01", "2032-06-30"),
    }

    if cycle not in CYCLE_DATES:
        raise ValueError(
            f"Cycle {cycle} not found. Available cycles: {list(CYCLE_DATES.keys())}"
        )
    return CYCLE_DATES[int(cycle)]
