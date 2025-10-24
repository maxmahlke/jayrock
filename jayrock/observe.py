import warnings

import json
import numpy as np
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation

import jayrock


class Observation:
    """JWST observation of a given target with a given instrument."""

    def __init__(self, target, instrument, config, report, date_obs, faint=False):
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
        faint : bool, optional
            Whether the observation was done at the faintest point. Default is False.
        """
        self.target = target
        self.instrument = instrument
        self.config = config
        self.report = report
        self.date_obs = date_obs
        self.faint = faint

        self.vmag = config["scene"][0]["spectrum"]["normalization"]["norm_flux"]
        self.thermal_flux = config["scene"][1]["spectrum"]["normalization"]["norm_flux"]

        self.wave, self.snr = report["1d"]["sn"]
        self.n_partial_saturated = report["1d"]["n_partial_saturated"][1]
        self.n_full_saturated = report["1d"]["n_full_saturated"][1]

        self.partial = np.array(self.n_partial_saturated) > 0
        self.full = np.array(self.n_full_saturated) > 0

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

    def plot_snr(self):
        """Plot the SNR as a function of wavelength."""
        jayrock.plotting.plot_snr(self)


def observe(target, instrument, faint=False):
    """Observe a target with a given instrument.

    Parameters
    ----------
    target : Target
        Target to observe.
    instrument : Instrument
        Instrument to use for the observation.
    faint : bool, optional
        Whether to observe asteroid at its faintest. Default is False.

    Returns
    -------
    Report
        Report of the observation.
    """
    jayrock.logging.logger.info(f"Observing {target} with {instrument}")

    if not hasattr(target, "visibility"):
        raise ValueError(
            "Target must have visibility information. Run target.query_visibility() first."
        )

    if not target.is_visible:
        raise ValueError(f"{target} is not visible. Stopping here.")

    # if not hasattr(target, "ephem"):
    #     raise ValueError(
    #         "Target must have ephemeris information. Run target.query_ephemerides() first."
    #     )

    if not hasattr(target, "thermal_flux"):
        raise ValueError(
            "Target must have thermal flux information. Run target.query_thermal_flux() first."
        )

    # ------
    # Build config
    config = build_default_calc(
        telescope="jwst", instrument=instrument.instrument, mode=instrument.mode
    )
    config["configuration"].update(instrument.config)

    # Set configuration elements decided by instrument
    if instrument.instrument == "nirspec":
        jayrock.logging.logger.info(
            f"[NIRSpec] Setting date_obs to when Vmag is {'minimum' if not faint else 'maximum'}.."
        )

        if not faint:
            date_obs = target.visibility.loc[
                target.visibility["V"] == target.vmag_min, "date_obs"
            ].values[0]
        else:
            date_obs = target.visibility.loc[
                target.visibility["V"] == target.vmag_max, "date_obs"
            ].values[0]

    elif instrument.instrument == "miri":
        jayrock.logging.logger.info(
            f"[MIRI] Setting date_obs to when thermal flux is {'minimum' if faint else 'maximum'}.."
        )

        if not faint:
            date_obs = target.visibility.date_obs.values[np.argmax(target.thermal_flux)]
        else:
            date_obs = target.visibility.date_obs.values[np.argmin(target.thermal_flux)]

    # Build scence from target properties
    scene = target.build_scene(date_obs)
    config["scene"] = scene

    jayrock.logging.logger.info(
        f"Observing on {date_obs} [V={target.vmag_date_obs:.2f}, Thermal: {target.thermal_flux_date_obs:.3f}mJy]"
    )

    # Strategy and background for SNR calculation
    config["strategy"]["reference_wavelength"] = instrument.reference_wavelength
    config["background_level"] = "none"

    # Sanity check
    # TODO: Insert link to docs
    if config["configuration"]["detector"]["ngroup"] < 5:
        jayrock.logging.logger.warning(
            f"Number of groups per integration is set to {config['configuration']['detector']['ngroup']}. "
            "It should be at least 5."
        )
    if config["configuration"]["detector"]["ngroup"] > 100:
        jayrock.logging.logger.warning(
            f"Number of groups per integration is set to {config['configuration']['detector']['ngroup']}. "
            "It should be at most 100. Increase the number of integrations per exposure instead."
        )

    # More info
    if instrument.instrument == "miri":
        grat = config["configuration"]["instrument"]["aperture"]
        disp = config["configuration"]["instrument"]["disperser"]
    elif instrument.instrument == "nirspec":
        grat = config["configuration"]["instrument"]["disperser"]
        disp = config["configuration"]["instrument"]["filter"]
    jayrock.logging.logger.info(
        f"{grat}|{disp} with {config['configuration']['detector']['ngroup']} groups per integration\n"
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = perform_calculation(config)

    return Observation(target, instrument, config, report, date_obs, faint=faint)


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
