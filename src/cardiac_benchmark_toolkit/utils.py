from dolfin import Mesh, MeshFunction, HDF5File, XDMFFile


def read_mesh(mesh_file: str) -> tuple[Mesh, MeshFunction]:
    tmp = mesh_file.split(".")  # [-1].split('.')
    file_type = tmp[-1]

    if file_type == "h5":
        mesh = Mesh()

        with HDF5File(mesh.mpi_comm(), mesh_file, "r") as hdf:
            hdf.read(mesh, "/mesh", False)
            boundaries = MeshFunction(
                "size_t", mesh, mesh.topology().dim() - 1
            )

            if hdf.has_dataset("boundaries"):
                hdf.read(boundaries, "/boundaries")
            else:
                if mesh.mpi_comm().Get_rank() == 0:
                    print(
                        "no <boundaries> datasets found in file {}".format(
                            mesh_file
                        )
                    )

    elif file_type == "xdmf":

        mesh = Mesh()

        with XDMFFile(mesh_file) as xf:
            xf.read(mesh)
            boundaries = MeshFunction(
                "size_t", mesh, mesh.topology().dim() - 1, 0
            )

            xf.read(boundaries)

    else:
        raise Exception("Mesh format not recognized. Use XDMF or HDF5.")

    return mesh, boundaries
