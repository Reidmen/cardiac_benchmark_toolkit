from importlib.metadata import PackageNotFoundError, version

__all__ = ["data", "ellipsoid_fiber_generation", "mesh_generation"]

try:
    __version__ = version("cardiac_benchmark_toolkit")
except PackageNotFoundError:
    __version__ = "UNDEFINED"
