import argparse
from pathlib import Path
import pathlib
import subprocess
import sys
import numpy as np

import meshio

from cardiac_benchmark_toolkit.data import DEFAULTS, MARKERS
from dolfin import (
    File,
    Function,
    FunctionSpace,
    HDF5File,
    Mesh,
    MeshFunction,
    MeshValueCollection,
    VectorFunctionSpace,
    XDMFFile,
)

from cardiac_benchmark_toolkit.ellipsoid_fiber_generation import read_mesh


def init_xdmf(
    path_to_folder: Path, filename: str = "fiber_directions"
) -> XDMFFile:
    """Initialize xdmf file using options dictionary"""
    path_to_file = path_to_folder.joinpath(f"{filename}.xdmf")
    xdmf_file = XDMFFile(str(path_to_file))
    xdmf_file.parameters["flush_output"] = True
    xdmf_file.parameters["functions_share_mesh"] = True

    return xdmf_file


def ellipsoid_mesh(
    hc: float,
    defaults: DEFAULTS = DEFAULTS(),
    markers: MARKERS = MARKERS(),
    path: str = "./meshes/",
):
    """Create truncated ellipsoid mesh with physiological dimensions.

    Surface IDs:
    1  endocardium
    2  epicardium
    3  base

    Args:
        hc: characteristic element size
        path: write path
    """
    GEO_CODE = f"""
        r_short_endo = {defaults.R_SHORT_ENDO};
        r_short_epi  = {defaults.R_SHORT_EPI};
        r_long_endo  = {defaults.R_LONG_ENDO};
        r_long_epi   = {defaults.R_LONG_EPI};
        quota_base = {defaults.QUOTA_BASE};

        mu_base = Acos(quota_base / r_long_endo);

        psize_ref = {hc};
        axisymmetric = 0;

        Geometry.CopyMeshingMethod = 1;
        Mesh.ElementOrder = 1;
        Mesh.Optimize = 1;
        Mesh.OptimizeNetgen = 1;
        Mesh.HighOrderOptimize = 1;

        Function EllipsoidPoint
            Point(id) = {{ r_long  * Cos(mu),
                           r_short * Sin(mu) * Cos(theta),
                           r_short * Sin(mu) * Sin(theta), psize }};
        Return

        center = newp; Point(center) = {{ 0.0, 0.0, 0.0 }};

        theta = 0.0;

        r_short = r_short_endo; r_long = r_long_endo;
        mu = -Pi;
        psize = psize_ref / 2.0;
        apex_endo = newp; id = apex_endo; Call EllipsoidPoint;
        mu = -1.0 * Acos(5.0 / 17.0);
        psize = psize_ref;
        base_endo = newp; id = base_endo; Call EllipsoidPoint;

        r_short = r_short_epi; r_long = r_long_epi;
        mu = -Pi;
        psize = psize_ref / 2.0;
        apex_epi = newp; id = apex_epi; Call EllipsoidPoint;
        mu = -1.0 * Acos(5.0 / 20.0);
        psize = psize_ref;
        base_epi = newp; id = base_epi; Call EllipsoidPoint;

        apex = newl; Line(apex) = {{ apex_endo, apex_epi }};
        base = newl; Line(base) = {{ base_endo, base_epi }};
        endo = newl;
        Ellipse(endo) = {{ apex_endo, center, apex_endo, base_endo }};
        epi  = newl;
        Ellipse(epi) = {{ apex_epi, center, apex_epi, base_epi }};

        ll1 = newll; Line Loop(ll1) = {{ apex, epi, -base, -endo }};
        s1 = news; Plane Surface(s1) = {{ ll1 }};

        If (axisymmetric == 0)
            sendoringlist[] = {{ }};
            sepiringlist[]  = {{ }};
            sendolist[] = {{ }};
            sepilist[]  = {{ }};
            sbaselist[] = {{ }};
            vlist[] = {{ }};

            sold = s1;
            For i In {{ 0 : 3 }}
                out[] = Extrude {{ {{ 1.0, 0.0, 0.0 }},
                                   {{ 0.0, 0.0, 0.0 }}, Pi/2 }}
                                {{ Surface{{sold}}; }};
                sendolist[i] = out[4];
                sepilist[i]  = out[2];
                sbaselist[i] = out[3];
                vlist[i] = out[1];
                bout[] = Boundary{{ Surface{{ sbaselist[i] }}; }};
                sendoringlist[i] = bout[1];
                sepiringlist[i] = bout[3];
                sold = out[0];
            EndFor

            // MYOCARDIUM
            Physical Volume(0) = {{ vlist[] }};
            // ENDO
            Physical Surface({markers.ENDOCARDIUM}) = {{ sendolist[] }};
            // EPI
            Physical Surface({markers.EPICARDIUM}) = {{ sepilist[] }};
            // BASE
            Physical Surface({markers.BASE}) = {{ sbaselist[] }};
            // ENDORING
            Physical Line(4) = {{ sendoringlist[] }};
            // EPIRING
            Physical Line(5) = {{ sepiringlist[] }};
        EndIf

        Physical Point("ENDOPT") = {{ apex_endo }};
        Physical Point("EPIPT") = {{ apex_epi }};
    """
    Path(path).mkdir(exist_ok=True, parents=True)

    geofile = Path(path).joinpath(f"ellipsoid_{hc}.geo")
    outfile = Path(path).joinpath(f"ellipsoid_{hc}.msh")
    with geofile.open("w") as f:
        f.write(GEO_CODE)
    subprocess.run(["gmsh", "-3", str(geofile)])

    # convert to dolfin XDMF
    mesh = meshio.read(str(outfile))

    xdmf_path = Path(path)
    xdmf_path.mkdir(exist_ok=True, parents=True)
    pth_tmp_msh = xdmf_path.joinpath("ellipsoid_meshio.xdmf")
    pth_tmp_bnd = xdmf_path.joinpath("ellipsoid_meshio_bound.xdmf")
    pth_aux_msh = xdmf_path.joinpath("ellipsoid_meshio.h5")
    pth_aux_bnd = xdmf_path.joinpath("ellipsoid_meshio_bound.h5")

    tetra_cells = mesh.get_cells_type("tetra")
    triangle_cells = mesh.get_cells_type("triangle")
    triangle_data = mesh.cell_data_dict["gmsh:physical"]["triangle"]

    tetra_mesh = meshio.Mesh(
        points=mesh.points, cells=[("tetra", tetra_cells)]
    )

    triangle_mesh = meshio.Mesh(
        points=mesh.points,
        cells=[("triangle", triangle_cells)],
        cell_data={"markers": [triangle_data]},
    )

    meshio.write(str(pth_tmp_msh), tetra_mesh, file_format="xdmf")

    meshio.write(str(pth_tmp_bnd), triangle_mesh, file_format="xdmf")

    mesh = Mesh()
    with XDMFFile(str(pth_tmp_msh)) as xf:
        xf.read(mesh)

    boundaries_mvc = MeshValueCollection("size_t", mesh)
    with XDMFFile(str(pth_tmp_bnd)) as xf:
        xf.read(boundaries_mvc)
    boundaries_mf = MeshFunction("size_t", mesh, boundaries_mvc)

    path_to_file = Path(path).joinpath(f"ellipsoid_{hc}.xdmf")

    with XDMFFile(str(path_to_file)) as xf:
        xf.write(mesh)
        xf.write(boundaries_mf)

    print("wrote {}".format(path))

    # remove temporary XDMF files and GEO/MSH
    geofile.unlink()
    # outfile.unlink()
    for pth_msh, pth_bnd in [
        (pth_tmp_msh, pth_tmp_bnd),
        (pth_aux_msh, pth_aux_bnd),
    ]:
        Path(pth_msh).unlink()
        Path(pth_bnd).unlink()


def biventricular_domain_from_mesh_and_fibers(
    path_to_folder: Path | str, element_space: int = 1
) -> None:
    """Process biventricular domain with with output fibers"""
    import vtk

    assert isinstance(element_space, int)
    assert isinstance(path_to_folder, (Path, str))
    print(f"Creating mesh with fibers in P{element_space}")

    xdmf_path = Path(path_to_folder)
    xdmf_path.mkdir(exist_ok=True, parents=True)

    domain_path = xdmf_path.joinpath("bi_ventricular.xdmf")
    fiber_path = xdmf_path.joinpath("fibers_biv.vtk")
    mesh, boundaries = read_mesh(domain_path.as_posix())

    V = VectorFunctionSpace(mesh, "CG", int(element_space))
    fiber = Function(V, name="fiber")
    sheet = Function(V, name="sheet")
    sheet_normal = Function(V, name="sheet_normal")

    (
        points_vtk,
        fiber_vtk,
        sheet_vtk,
        sheet_normal_vtk,
    ) = load_vtk_data_from_filepath(fiber_path.as_posix())

    names = ["fiber", "sheet", "sheet_normal"]
    directions = [fiber, sheet, sheet_normal]
    vtk_directions = [fiber_vtk, sheet_vtk, sheet_normal_vtk]
    idx_0, _, _ = compute_transfer_indices_vtk_to_dolfin(V, points_vtk)

    load_fibers_with_vtk_data(directions, vtk_directions, names, idx_0)
    save_mesh_and_boundaries_into_vtk(xdmf_path, mesh, boundaries)
    save_fibers_into_h5_and_vtk(xdmf_path, mesh, directions, names)
    save_into_xdmf(xdmf_path, fiber, sheet, sheet_normal)


def save_fibers_into_h5_and_vtk(
    pathfolder: pathlib.Path,
    mesh: Mesh,
    directions: list[Function],
    names: list[str],
) -> None:
    print("Saving fiber files in the standard format")
    for func, name in zip(directions, names):
        hdf5_file = pathfolder.joinpath(
            f"fibers/h5_format/bi_ventricular_{name}.h5"
        )

        with HDF5File(mesh.mpi_comm(), hdf5_file.as_posix(), "w") as hdf:
            hdf.write(func, "/" + name)

        vtk_fiber = File(
            pathfolder.joinpath(
                f"fibers/pvd_format/bi_ventricular_{name}.pvd"
            ).as_posix()
        )
        vtk_fiber << func
        print(
            f"Wrote {name} in path {hdf5_file.as_posix()} alongside VTK format"
        )


def load_fibers_with_vtk_data(
    directions: list[Function],
    vtk_directions: list[np.ndarray],
    names: list[str],
    idx_map: list[int],
) -> None:
    print("Loading fiber files in the standard format")
    for func, vtk_func, name in zip(directions, vtk_directions, names):
        dofs_0 = func.function_space().sub(0).dofmap().dofs()
        dofs_1 = func.function_space().sub(1).dofmap().dofs()
        dofs_2 = func.function_space().sub(2).dofmap().dofs()

        func.vector()[dofs_0] = vtk_func[:, 0][idx_map]
        func.vector()[dofs_1] = vtk_func[:, 1][idx_map]
        func.vector()[dofs_2] = vtk_func[:, 2][idx_map]


def save_mesh_and_boundaries_into_vtk(
    folderpath: pathlib.Path, mesh: Mesh, boundaries: MeshFunction
) -> None:
    assert isinstance(mesh, Mesh)
    # boundaries type MeshFunctionSizet (dolfin.cpp)
    vtkfile = File(folderpath.joinpath("vtk/bi_ventricular.pvd").as_posix())
    mfunc = MeshFunction("size_t", mesh, 3)
    mfunc.set_all(0)

    vtkfile << mesh
    vtkfile << boundaries
    vtkfile << mfunc

    print("Wrote mesh, boundaries in PVD format")


def save_into_xdmf(
    path_to_save: Path,
    fiber: Function,
    sheet: Function,
    sheet_normal: Function,
) -> None:
    print(f"Saving into XDMF format in {path_to_save.as_posix()}")
    xdmf_file = init_xdmf(path_to_save)
    xdmf_file.write(fiber, 0.0)
    xdmf_file.write(sheet, 0.0)
    xdmf_file.write(sheet_normal, 0.0)
    xdmf_file.close()


def load_vtk_data_from_filepath(
    filepath: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Loads vtk points, fiber, sheet and sheet-normal directions"""
    import vtk

    print("Loading vtk files")
    if isinstance(filepath, str):
        filepath = Path(filepath)

    fibers_fields = vtk.vtkGenericDataObjectReader()
    fibers_fields.SetFileName(filepath.as_posix())
    fibers_fields.ReadAllVectorsOn()
    fibers_fields.Update()

    points = np.array(fibers_fields.GetOutput().GetPoints().GetData())
    fiber_vtk = np.array(
        fibers_fields.GetOutput().GetPointData().GetArray("f0")
    )
    sheet_vtk = np.array(
        fibers_fields.GetOutput().GetPointData().GetArray("s0")
    )
    sheet_normal_vtk = np.array(
        fibers_fields.GetOutput().GetPointData().GetArray("n0")
    )

    return points, fiber_vtk, sheet_vtk, sheet_normal_vtk


def compute_transfer_indices_vtk_to_dolfin(
    V: FunctionSpace, vtk_points: np.ndarray
) -> tuple[list[int], list[int], list[int]]:
    """Computers transfer indices from vtk point arrays to dolfin function space"""
    print("Computing transfer index between vtk and dolfin")
    coords_0, coords_1, coords_2 = get_subspace_coordinates(V)
    idx_0, idx_1, idx_2 = [], [], []

    for coords, transfer_idx in zip(
        [coords_0],
        [idx_0],
    ):
        for point in coords:
            distance_vector = np.sqrt(
                np.sum((point - vtk_points) ** 2, axis=1)
            )
            min_index = np.argmin(distance_vector)
            transfer_idx.append(min_index)

    return idx_0, idx_0, idx_0


def get_subspace_coordinates(
    V: FunctionSpace,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Get coordinates of subspaces (dim 3 fixed)"""
    total_dofs_coordinates = V.tabulate_dof_coordinates()
    dofs_0 = V.sub(0).dofmap().dofs()
    dofs_1 = V.sub(1).dofmap().dofs()
    dofs_2 = V.sub(2).dofmap().dofs()

    coords_sub_0 = total_dofs_coordinates[dofs_0, :]
    coords_sub_1 = total_dofs_coordinates[dofs_1, :]
    coords_sub_2 = total_dofs_coordinates[dofs_2, :]

    return coords_sub_0, coords_sub_1, coords_sub_2


def get_parser():
    """Get arguments parser.

    Returns:
        parser
    """
    parser = argparse.ArgumentParser(
        description="""
        Generate ellipsoid mesh, with optional path and element size.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-create_from_data",
        "--create_mesh_and_fibers_from_folder",
        type=str,
        default=None,
        help="Path to mesh and fiber (single folder)",
    )
    parser.add_argument(
        "-deg",
        "--element_degree",
        type=int,
        default=0,
        help="FE degree for interpolation of fibers",
    )
    parser.add_argument(
        "-path",
        "--path_to_save",
        type=str,
        default="./meshes/",
        help="Path where meshes are written",
    )
    parser.add_argument(
        "-size",
        "--element_size",
        metavar="hc",
        type=float,
        default=0.005,
        help="Truncated ellipsoid mesh with characteristic mesh size",
    )
    return parser


if __name__ == "__main__":
    args: argparse.Namespace = get_parser().parse_args()
    print(args)

    if len(sys.argv) > 1:
        if (
            args.create_mesh_and_fibers_from_folder is not None
            and args.element_degree > 0
        ):
            biventricular_domain_from_mesh_and_fibers(
                args.create_mesh_and_fibers_from_folder, args.element_degree
            )
        else:
            ellipsoid_mesh(args.element_size, path=args.path_to_save)
