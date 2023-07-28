import subprocess
import meshio
from pathlib import Path
import argparse
import sys
from dolfin import XDMFFile, Mesh, MeshFunction, MeshValueCollection

from cardiac_benchmark_toolkit.data import MARKERS, DEFAULTS


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

    if len(sys.argv) > 1:
        ellipsoid_mesh(args.element_size, path=args.path_to_save)
