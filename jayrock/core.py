from astropy.time import Time
from pandeia.engine.calc_utils import build_default_calc
import rocks
from rich import print

from . import jwst
from . import target_model
import jayrock


class Observation:
    """JWST observation of a given target with a given instrument."""

    def __init__(
        self,
        target,
        instrument,
        mode,
        cycle=None,
        date_start=None,
        date_end=None,
        config=None,
        label=None,
    ):
        """Configure the observation.

        Parameters
        ----------
        target : str
            Name of the target to observe. Passed to rocks.id for identification.
        instrument : str
            Name of the instrument to use for the observation.
        mode : str
            Mode of the instrument to use for the observation.
        cycle : int, optional
            JWST observation cycle to determine start and end dates of observation.
            Must be provided if date_start and date_end are not given. Default is None.
        date_start : str, optional
            Start date of the observation in ISO format (YYYY-MM-DD). Can be
            omitted if cycle is given. Default is None.
        date_end : str, optional
            End date of the observation in ISO format (YYYY-MM-DD). Can be
            omitted if cycle is given. Default is None.
        config : dict, optional
            Configuration to pass to the JWST ETC calculation. See

            https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial/pandeia-quickstart
        """

        # Set target
        self.target = rocks.Rock(target)

        # Set dates
        self.cycle = cycle

        if self.cycle is not None:
            self.date_start, self.date_end = jwst.get_cycle_dates(cycle)
        elif date_start is not None and date_end is not None:
            self.date_start, self.date_end = date_start, date_end
        else:
            raise ValueError(
                "Either cycle or both date_start and date_end must be provided."
            )

        self.date_start = Time(self.date_start, format="iso", scale="utc")
        self.date_end = Time(self.date_end, format="iso", scale="utc")

        # Set instrument configuration
        self.instrument = instrument.lower()
        self.mode = mode.lower()

        # self.config = build_default_calc("jwst", self.instrument, self.mode)
        self.config = config if config is not None else {}
        self.file_suffix = f"_{label}" if label is not None else ""
        self.label = label

        self.create_output_dirs()

    def create_output_dirs(self):
        """Create output directories if they don't exist."""
        self.PATH_TARGET = jayrock.PATH_OUT / self.target.filename
        self.PATH_TARGET.mkdir(parents=True, exist_ok=True)

        self.PATH_INST = self.PATH_TARGET / self.instrument
        self.PATH_INST.mkdir(parents=True, exist_ok=True)

    def store_visibility(self):
        """Store visibility windows to a CSV file."""
        PATH_VIS = self.PATH_TARGET / f"vis.csv"
        self.vis.to_csv(PATH_VIS, index=False)
