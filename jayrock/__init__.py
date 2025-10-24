"""For JWST observations of space rocks."""

# Welcome to jayrock
__version__ = "0.1"

from . import neatm, plotting
from .logging import set_log_level  # noqa
from .instrument import Instrument
from .observe import get_cycle_dates, Observation, observe
from .target import Target
