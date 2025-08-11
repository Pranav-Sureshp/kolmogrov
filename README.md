# Kolmogrov length scale and dissipation

Add Openfoam code to calculate Kolmogrov scale and dissipation

# Purpose

This code is a basic OpenFOAM utility that reads two fields from a case directory:

    Velocity field U (a volVectorField)
    Turbulent viscosity nut (a volScalarField)

It sets up the simulation environment and reads these fields from disk using OpenFOAM's object-oriented framework.

Header Files Explained
#include "fvCFD.H"

This is a master header for finite volume CFD applications. It includes:

    Mesh handling (fvMesh)
    Field types (volScalarField, volVectorField, etc.)
    IO operations (IOobject)
    Time control (Time)
    Logging (Info)
    Mathematical operations and solvers

This header is typically used in solver and utility development.
#include "setRootCase.H"

Sets the root case directory based on command-line arguments. It initializes the OpenFOAM environment and ensures the correct case folder is used.
#include "createTime.H"

Creates the Time object, which manages simulation time and controls reading/writing of time directories. This is essential for accessing time-dependent data.
#include "createMesh.H"

Creates the fvMesh object, which represents the computational mesh. It reads mesh data from the constant/polyMesh directory and initializes the mesh for field operations.


- Reading Fields
- Velocity Field U
    
    volVectorField: A volume field of vectors (used for velocity).
    IOobject constructor:
  -      "U": Name of the field file to read (U).
  -      runTime.timeName(): Current time directory (e.g., "0", "100").
  -      mesh: Mesh associated with the field.
  -      IOobject::MUST_READ: Field must be read from disk.
  -      IOobject::AUTO_WRITE: Field will be automatically written during output.

This reads the velocity field from the current time directory.

# Gradient and Derived Fields
# Purpose

This code computes:

    The velocity gradient tensor (gradU)
    The symmetric part of the velocity gradient (S)
    The effective dissipation rate (epsilon_eff)
    The Kolmogorov velocity scale (kolmogorov)

These quantities are useful in turbulence modeling, flow diagnostics, and post-processing.

-    fvc::grad(U): Computes the gradient of the velocity field U using finite volume calculus (explicit method).
-    volTensorField: A volume field of tensors, representing the spatial derivatives of the velocity vector field.

-    symm(gradU): Extracts the symmetric part of the tensor field gradU, i.e., $S=0.5(\nabla U + \nabla U^T)$
-    volSymmTensorField: A volume field of symmetric tensors, often used to compute strain rate tensors in turbulence models.
- Dissipation : $2(\nu + \nu_t) S_{ij} S_{ij}$.
- Kolmogrov length scale : $\eta = (\nu^3/\epsilon)^{1/4} $.

