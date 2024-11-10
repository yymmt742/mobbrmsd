from .dataclass import molecules, molecular_system
from .mobbrmsd import mobbrmsd
from .demo._coord_generator import coord_generator
from importlib.metadata import version

__version__ = version(__package__)
