![MIT](https://img.shields.io/badge/License-MIT-green)
![Black](https://img.shields.io/badge/Style-Black-black)
# Cardiac Benchmark Toolkit

This repository contains basic scripts that allow you to reproduce the mesh
as well as the fibers in the cardiac mechanics benchmark in dolfin (legacy). 

*For the [dolfinx](https://github.com/FEniCS/dolfinx) implementation of the tools, [check here](https://github.com/Reidmen/cardiac_benchmark_toolkitx)*

## Installation

*Docker*
Run the following command to start a container with all the required dependencies:

```shell
docker run --name dolfin-stable -v $(pwd):/home/shared -w /home/shared -ti ghcr.io/scientificcomputing/fenics-gmsh:2023-04-21
```

In order to enter the shell, use:

```shell
docker exec -ti dolfin-stable /bin/bash -l
```

**Note** docker image is cortesy of *Simula Lab*.

## Quickstart

### Fiber Generation
This repository provides an ellispoid tagged mesh for reference purposes in the folder `./meshes`.
Use the mesh `./meshes/ellipsoid_0.005.xdmf`, you can create the fibers as follows:

```shell
cardiac_benchmark_toolkit/ellipsoid_fiber_generation.py ./meshes/ellipsoid.xdmf
```

If succesfull, the script will create fibers in `xdmf` and `vtk` files in a `./results/` folder.

Further options can be found with:

```shell
cardiac_benchmark_toolkit/ellipsoid_fiber_generation.py --help
```


## Mesh Generation
The ellipsoid domain can be created with a characteristic element size (default is `0.005 [m]`) using the
script `mesh_generation.py`. It will create a folder structure `./results` to stores `xdmf` as well as `pvd` formats
for further usage.

To execute it, consider the following command:
```shell
cardiac_benchmark_toolkit/mesh_generation.py -size 0.007
```

It will create an ellipsoid mesh with characteristic element size of `0.007 [m]`. You can use it in conjuntion with the
`ellipsoid_fiber_generation` to create the fiber directions for your specific simulation (and benchmark).
