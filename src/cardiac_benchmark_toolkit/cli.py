from importlib.metadata import metadata
import pathlib

import rich_click as click

from cardiac_benchmark_toolkit.ellipsoid_fiber_generation import (
    create_fibers_for_ellipsoid,
)
from cardiac_benchmark_toolkit.mesh_generation import (
    create_biventricular_domain_from_mesh_and_fibers,
    create_ellipsoid_mesh,
)

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
@click.argument("outdir", required=True, type=click.Path())
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
    create_ellipsoid_mesh(char_length, path=str(outdir))


@click.command(help="Create biv ellipsoid and fibers from data")
@click.argument("datadir", required=True, type=click.Path())
@click.option(
    "--element_degree",
    type=int,
    default=1,
    help="Finite Element degree for interpolation of fibers",
)
def create_biv_domain_from_mesh_and_fibers(
    datadir: pathlib.Path, element_degree: int
) -> None:
    create_biventricular_domain_from_mesh_and_fibers(
        str(datadir), element_degree
    )


@click.command(
    help="Create fiber, sheet and sheet_normal directions for ellipsoid"
)
@click.argument("path_to_mesh", required=True, type=click.Path())
@click.option(
    "--element_degree",
    type=int,
    default=1,
    help="Finite Element degree for interpolation of fibers",
)
@click.option(
    "--outdir",
    type=click.Path(),
    default=pathlib.Path("./results"),
    help="Directory to output data, default is './results/'",
)
def create_fibers_for_ellipsoid_mesh(
    path_to_mesh: pathlib.Path, element_degree: int, outdir: pathlib.Path
) -> None:
    if isinstance(path_to_mesh, str):
        path_to_mesh = pathlib.Path(path_to_mesh)

    if isinstance(outdir, str):
        outdir = pathlib.Path(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    assert isinstance(path_to_mesh, pathlib.Path)
    assert isinstance(outdir, pathlib.Path)
    assert path_to_mesh.is_file()
    assert outdir.is_dir()

    create_fibers_for_ellipsoid(
        path_to_mesh.as_posix(), element_degree, outdir.as_posix()
    )


app.add_command(create_ellipsoid)
app.add_command(create_biv_domain_from_mesh_and_fibers)
app.add_command(create_fibers_for_ellipsoid_mesh)
