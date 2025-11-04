import json

import numpy as np
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.instrument_factory import InstrumentFactory

import jayrock

URLS = {
    "nirspec": {
        "detector": "https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-observing-strategies/nirspec-detector-recommended-strategies",
        "dither": "https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-operations/nirspec-dithers-and-nods",
    },
    "miri": {
        "dither": "https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-operations/miri-dithering"
        "detector": "https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-strategies"
    },
}


class Detector:
    """Detector of observation."""

    def __init__(self, config, instrument):
        """
        Parameters
        ----------
        config : dict
            instrument['detector'] pandeia configuration dictionary.
        instrument : Instrument
            Instrument instance that this detector belongs to.
        """
        self._ngroup = config["ngroup"]
        self._nint = config["nint"]
        self._nexp = config["nexp"]
        self._readout_pattern = config["readout_pattern"]
        self._subarray = config["subarray"]

        self.instrument = instrument

    @property
    def ngroup(self):
        return self._ngroup

    @ngroup.setter
    def ngroup(self, ngroup):
        if ngroup < 5:
            jayrock.logging.logger.warning(
                f"Setting 'ngroup={ngroup}'. 'ngroup' values below 5 are discouraged."
                f" More information:\n{URLS[self.instrument.instrument]['detector']}"
            )
        if ngroup > 100 and self.instrument.instrument == "nirspec":
            jayrock.logging.logger.warning(
                f"Setting 'ngroup={ngroup}'. 'ngroup' values larger than 100 are discouraged. Consider increasing 'nint' instead."
                f" More information:\n{URLS[self.instrument.instrument]['detector']}"
            )
        if self.readout_pattern == "nrsirs2rapid" and self.nint * ngroup > 1024:
            jayrock.logging.logger.warning(
                'nint * ngroup has to be <= 1024 for "nrsirs2rapid" readout pattern.'
            )
        if isinstance(ngroup, np.int64):
            ngroup = int(ngroup)
        self._ngroup = ngroup

    @property
    def nint(self):
        return self._nint

    @nint.setter
    def nint(self, nint):
        if self.readout_pattern == "nrsirs2rapid" and nint * self.ngroup > 1024:
            jayrock.logging.logger.warning(
                'nint * ngroup has to be <= 1024 for "nrsirs2rapid" readout pattern.'
            )
        if isinstance(nint, np.int64):
            nint = int(nint)
        self._nint = nint

    @property
    def subarray(self):
        return self._subarray

    @subarray.setter
    def subarray(self, subarray):
        self._subarray = subarray

    @property
    def nexp(self):
        return self._nexp

    @nexp.setter
    def nexp(self, nexp):
        if nexp == 1:
            jayrock.logging.logger.warning(
                f"Setting 'nexp={nexp}'. Dithering ('nexp' > 1) is generally recommended."
                f" More information:\n{URLS[self.instrument.instrument]['dither']}"
            )

        self._nexp = nexp

    @property
    def readout_pattern(self):
        return self._readout_pattern

    @readout_pattern.setter
    def readout_pattern(self, rdpt):
        rdpt = rdpt.lower()

        RDPT_NIRSPEC = ["nrs", "nrsrapid", "nrsirs2", "nrsirs2rapid"]

        if self.instrument.instrument == "nirspec":
            if rdpt not in RDPT_NIRSPEC:
                raise ValueError(
                    f"Invalid 'readout_pattern' for NIRSpec: {rdpt}. Choose from {RDPT_NIRSPEC}."
                )

            if rdpt != "nrsirs2rapid":
                jayrock.logging.logger.warning(
                    f"Setting 'readout_pattern={rdpt}'. The 'nrsirs2rapid' readout pattern is "
                    "recommended for NIRSpec observations."
                    f" More information:\n{URLS[self.instrument.instrument]['detector']}"
                )

        self._readout_pattern = rdpt


class Instrument:
    """Instrument of observation."""

    def __init__(self, instrument, mode, config=None):
        """
        Parameters
        ----------
        instrument : str
            Instrument name (e.g., 'nirspec', 'miri').
        mode : str
            Observation mode (e.g., 'ifu', 'mrs').
        config : dict, optional
            Configuration dictionary for the instrument. If None, a default configuration will be used.
            Configuration dictionary should be the 'configuration' level of the pandeia calculation config.
        """
        self._instrument = instrument.lower()
        self._mode = mode.lower()

        if config is not None:
            self._config = config
        else:
            self._config = build_default_calc(
                telescope="jwst", instrument=self._instrument, mode=self._mode
            )["configuration"]

        if self._instrument == "nirspec":
            # set better default
            self._config["instrument"]["readout_pattern"] = "nrsirs2rapid"

        self.detector = Detector(self._config["detector"], self)

    def __repr__(self):
        return f"Instrument(instrument={self.instrument}, mode={self.mode})"

    @property
    def config(self):
        """Pandeia instrument config dictionary."""
        for attr in ["aperture", "disperser", "filter"]:
            self._config["instrument"][attr] = getattr(self, attr)
        for attr in ["ngroup", "nint", "nexp", "readout_pattern", "subarray"]:
            self._config["detector"][attr] = getattr(self.detector, attr)
        return self._config

    @property
    def instrument(self):
        """Write-protected instrument setting."""
        return self._instrument

    @property
    def mode(self):
        """Write-protected mode setting."""
        return self._mode

    @property
    def disperser(self):
        """Disperser used in instrument."""
        return self._config["instrument"]["disperser"]

    @disperser.setter
    def disperser(self, disp):
        self._config["instrument"]["disperser"] = disp.lower()

    @property
    def filter(self):
        """Filter used in instrument."""
        return self._config["instrument"]["filter"]

    @filter.setter
    def filter(self, filt):
        self._config["instrument"]["filter"] = filt.lower()

    @property
    def aperture(self):
        """Aperture used in instrument."""
        return self._config["instrument"]["aperture"]

    @aperture.setter
    def aperture(self, apt):
        self._config["instrument"]["aperture"] = apt.lower()

    @property
    def name(self):
        return {"nirspec": "NIRSpec", "miri": "MIRI"}[self._instrument]

    @property
    def texp(self):
        """Total exposure time in seconds."""
        inst = InstrumentFactory(config=self.config)
        return inst.the_detector.exposure_spec.total_exposure_time

    @property
    def primary(self):
        """Return primary optical element of the instrument."""
        return {"nirspec": self.disperser, "miri": self.aperture}[self.instrument]

    @property
    def secondary(self):
        """Return secondary optical element of the instrument."""
        return {"nirspec": self.filter, "miri": self.disperser}[self.instrument]

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
            }[self.filter]

        if self.instrument == "miri":
            if self.mode == "mrs":
                return {
                    "ch1": {"short": 5.335, "medium": 6.16, "long": 7.109},
                    "ch2": {"short": 8.154, "medium": 9.415, "long": 10.85},
                    "ch3": {"short": 12.504, "medium": 14.5, "long": 16.745},
                    "ch4": {"short": 19.29, "medium": 22.485, "long": 26.22},
                }[self.aperture][self.disperser]
            elif self.mode in ["lrsslit", 'lrsslitless']:
                return 10.7

    def print_config(self):
        """Print the observation configuration."""
        print(json.dumps(self.config, indent=4))


    def set_snr_target(self, snr, target, date_obs, wave=None, bounds=None):
        """Configure detector settings to reach a target SNR.

        Parameters
        ----------
        snr : float
            The SNR to reach.
        target : jayrock.Target
            The observation target.
        date_obs : str
            The observation date string in YYYY-MM-DD format.
        wave : float, optional
            Wavelength in microns at which to optimise the SNR. If None, the
            reference wavelength of the instrument is used.
        bounds : dict, optional
            Optional bounds for detector settings. Valid keys are 'ngroup' and 'nint'.
            Each key should map to a tuple of (min, max) values.

        Returns
        -------
        bool
            True if the target SNR or higher was achieved, False otherwise (e.g. saturation).

        Notes
        -----
        This method modifies the detector's 'ngroup' and 'nint' settings to
        achieve the desired SNR using a binary search algorithm. It prioritises
        minimising 'nint'.
        """
        # Determine if we are searching for a single target or a range
        is_range_target = isinstance(snr, (list, tuple))
        if is_range_target:
            min_snr, max_snr = snr

        wave = wave if wave is not None else self.reference_wavelength
        config_override = {"strategy": {"reference_wavelength": wave}}
        scene = target.build_scene(date_obs)

        # Define search space with defaults
        bounds = bounds if bounds is not None else {}
        ngroup_min, ngroup_max = bounds.get("ngroup", (5, 100))
        nint_min, nint_max = bounds.get("nint", (1, 100))

        # TODO: Separate the logic here into SNR range vs single SNR target functions

        if is_range_target:
            jayrock.logging.logger.info(
                f"Searching for minimum ngroup|nint to reach SNR range {min_snr:.1f}-{max_snr:.1f} at {wave:.2f}μm"
            )
        else:
            jayrock.logging.logger.info(
                f"Searching for minimum ngroup|nint to reach SNR {snr:.1f} at {wave:.2f}μm"
            )

        def run_obs(nint: int, ngroup: int):
            """Helper function to run an observation and check for saturation."""
            self.detector.nint = nint
            self.detector.ngroup = ngroup

            obs = jayrock.observe(
                target,
                self,
                date_obs=date_obs,
                scene=scene,
                config_override=config_override,
                verbose=False,
            )

            if any(obs.fully_saturated):
                jayrock.logging.logger.warning(
                    f"  nint={nint} | ngroup={ngroup:<3} -> FULL SATURATION."
                )
                # Return None to signal saturation
                return None

            jayrock.logging.logger.info(
                f"  nint={nint} | ngroup={ngroup:<3} -> SNR={obs.snr:<5.1f} | Texp={obs.texp / 60:.1f}min"
            )
            return obs

        # ------
        # Run search
        for nint in range(nint_min, nint_max + 1):
            # Reset ngroup bounds for the binary search
            ngroup_low = ngroup_min
            ngroup_high = ngroup_max

            # ------
            # Check 1: Does the minimum ngroup saturate? Does it reach the target SNR?
            if nint == nint_min:
                obs = run_obs(nint, ngroup_min)

                if obs is None:
                    jayrock.logging.logger.error(
                        f"Minimum nint/ngroup saturated. Stopping search. Providing parameter 'bounds' might avoid this."
                    )
                    return False

                if is_range_target:
                    if obs.snr >= min_snr:
                        # Min ngroup already sufficient, done
                        self.detector.nint = nint
                        self.detector.ngroup = ngroup_min

                        jayrock.logging.logger.info(
                            f"Done. Setting ngroup={ngroup_min}, nint={nint}."
                        )
                        self.estimated_snr = obs.snr
                        return True
                else:
                    if obs.snr >= snr:
                        # Min ngroup already sufficient, done
                        self.detector.nint = nint
                        self.detector.ngroup = ngroup_min

                        jayrock.logging.logger.info(
                            f"Done. Setting ngroup={ngroup_min}, nint={nint}."
                        )
                        self.estimated_snr = obs.snr
                        return True

            # ------
            # Check 2: Can we reach the target SNR with this nint?
            obs = run_obs(nint, ngroup_max)

            if (
                obs is not None
            ):  # if not saturated. if it is, we might have more luck with lower ngroup
                if is_range_target:
                    if min_snr <= obs.snr <= max_snr:
                        # Max ngroup sufficient, done
                        self.detector.nint = nint
                        self.detector.ngroup = ngroup_max

                        jayrock.logging.logger.info(
                            f"Done. Setting ngroup={ngroup_max}, nint={nint}."
                        )
                        self.estimated_snr = obs.snr
                        return True
                    if obs.snr < min_snr:
                        continue  # Max ngroup is still too low, try next nint
                else:
                    if obs.snr < snr:
                        continue  # Max ngroup is still too low, try next nint

            # Binary Search on ngroup: Find the minimum ngroup that yields SNR >= snr
            while ngroup_low < ngroup_high:
                # Check for convergence protection
                if ngroup_low == ngroup_high - 1:
                    # If only two discrete values remain, test the lower one first.
                    # This ensures we find the *minimum* that satisfies the condition.
                    obs = run_obs(nint, ngroup_low)

                    if obs is None:
                        ngroup_high = ngroup_low
                        break  # use the lower value
                    if is_range_target:
                        if min_snr <= obs.snr <= max_snr:
                            # Max ngroup sufficient, done
                            self.detector.nint = nint
                            self.detector.ngroup = ngroup_low
                            self.estimated_snr = obs.snr

                            jayrock.logging.logger.info(
                                f"Done. Setting ngroup={ngroup_low}, nint={nint}."
                            )
                        return True
                        if obs.snr >= max_snr:
                            ngroup_high = ngroup_low
                            break  # use the lower value
                    else:
                        if obs.snr >= snr:
                            ngroup_high = ngroup_low
                            break  # use the lower value
                    break  # use the higher value

                # Try mid-point ngroup
                ngroup_mid = (ngroup_low + ngroup_high) // 2
                obs = run_obs(nint, ngroup_mid)

                if is_range_target:
                    if min_snr <= obs.snr <= max_snr:
                        # Max ngroup sufficient, done
                        self.detector.nint = nint
                        self.detector.ngroup = ngroup_mid
                        self.estimated_snr = obs.snr

                        jayrock.logging.logger.info(
                            f"Done. Setting ngroup={ngroup_mid}, nint={nint}."
                        )
                        return True
                    if obs.snr < min_snr:
                        # SNR too low, increase ngroup. Search upper half of
                        # current interval
                        ngroup_low = ngroup_mid + 1
                    else:
                        # SNR is sufficient, try lower ngroup to minimise time.
                        # Search lower half of current interval
                        ngroup_high = ngroup_mid
                else:
                    if obs.snr < snr:
                        # SNR too low, increase ngroup. Search upper half of
                        # current interval
                        ngroup_low = ngroup_mid + 1
                    else:
                        # SNR is sufficient, try lower ngroup to minimise time.
                        # Search lower half of current interval
                        ngroup_high = ngroup_mid

            # Done. Update detector settings
            self.detector.nint = nint
            self.detector.ngroup = ngroup_high

            jayrock.logging.logger.info(
                f"Done. Setting ngroup={ngroup_high}, nint={nint}."
            )
            self.estimated_snr = obs.snr
            return True
