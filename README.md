# 3d_BEM
The 3d_BEM program is a tool for solving the Helmholtz equation in acoustics, using the Boundary Element Method. The BEM solver uses H-matrices and adaptive cross approximation to speed up the computation. Acoustic simulations can be controlled with a simple script interpreter. Mesh files (GMSH-format) can be supplied to set up complex geometries with general Robin boundary conditions (impedance, radiation, reflection, etc.). The geometries and the sound field solutions are visualized with a 3d viewer.

## Compiling 3d_BEM
Some instructions on how to compile 3d_BEM from source code.

###### Requirements
Before you can compile 3d_BEM, you need to have the following software installed on your computer:

 - A C++20-compliant compiler, such as GCC 8 or higher.
 - The Qt library, version 6.0 or higher .
 - Eigen 3.3.9 or higher, for linear algebra operations (https://eigen.tuxfamily.org). Can be installed under Ubuntu with the following command:
```
 sudo apt-get install libeigen3-dev
```
 - OpenGL libraries -lGL -lGLU; can be installed under Ubuntu with the following command:

```
 sudo apt-get install libgl1-mesa-dev libglu1-mesa-dev 
```
###### Download
The first step is to download the 3d_BEM source code from the project's website or from its Git repository.

###### Building
To compile 3d_BEM, follow these steps:

 1. Open a command prompt or terminal window and navigate to the desired build directory.
 2. Call qmake, referencing the .pro file in the source directory.
 3. Call make.

```
 * cd <build directory> 
 * qmake <source directory>
 * make
```
Alternatively just open the .pro file with the QtCreator and compile from there.

## How to use 3d_BEM
Some example 3d_BEM projects are provided with the source code in the *Examples* folder.
To be continued...

## Literature
Exhaustive background information on the Boundary Element Method and its acceleration via H-matrices can be found in the *Literature* folder.
