from importlib.metadata import PackageNotFoundError, version
from . import (
    data,
    ellipsoid_fiber_generation,
    mesh_generation,
)

__all__ = ["data", "ellipsoid_fiber_generation", "mesh_generation"]

try:
    __version__ = version("cardiac_benchmark_toolkit")
except PackageNotFoundError:
    __version__ = "UNDEFINED"
