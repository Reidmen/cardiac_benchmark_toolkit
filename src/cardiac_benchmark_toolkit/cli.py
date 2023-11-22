from importlib.metadata import metadata
import pathlib

from mesh_generation import (
    ellipsoid_mesh,
    biventricular_domain_from_mesh_and_fibers,
)

import rich_click as click

meta = metadata("cardiac_benchmark_toolkit")
__version__ = meta["version"]
__author__ = meta["author"]
__licence__ = meta["license"]


@click.group()
@click.version_option(__version__, prog_name="cardiac_benchmark_toolkit")
def app() -> None:
    """Cardiac Benchmark Toolkit - Basic set of tools for benchmarking"""
    pass


@click.command(help="Create ellipsoid geometry")
@click.argument(
    "outdir", required=True, type=click.Path(), help="Path to save mesh"
)
@click.option(
    "--char-length",
    default=0.01,
    type=float,
    help="Characteristic length of mesh",
    show_default=True,
)
def create_ellipsoid(
    outdir: pathlib.Path,
    char_length: float = 0.01,
) -> None:
    ellipsoid_mesh(char_length, path=str(outdir))

@click.command(help="Create biv ellipsoid and fibers from data")
@click.argument(
    "datadir", required=True, type=click.Path(), help="directory with data"
)
@click.option(
    "--element_degree",
    type=int,
    default=1,
    help="Finite Element degree for interpolation of fibers",
)
def create_biv_domain_from_mesh_and_fibers(
    datadir: pathlib.Path, element_degree: int
) -> None:
    biventricular_domain_from_mesh_and_fibers(str(datadir), element_degree)


app.add_command(create_ellipsoid)
app.add_command(create_biv_domain_from_mesh_and_fibers)
