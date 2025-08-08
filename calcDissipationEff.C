#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Reading velocity field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading turbulent viscosity nut\n" << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Set kinematic molecular viscosity (update to match your case!)
    const scalar nu = 3.4e-7;

    Info << "Calculating velocity gradient tensor (gradU)\n" << endl;
    volTensorField gradU = fvc::grad(U);
    volSymmTensorField S = symm(gradU);

    Info << "Computing effective dissipation rate (epsilon_eff) cell-by-cell..." << endl;
    volScalarField epsilon_eff
    (
        IOobject
        (
            "epsilon_eff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    );

    volScalarField kolmogorov
    (
        IOobject
        (
            "kolmogorov",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,1,0,0,0,0,0), 0.0)
    );

    // Main for-loop: Loop through all cells and assign epsilon_eff per cell
    forAll(mesh.cells(), celli)
    {
        scalar nut_cell = nut[celli];       // turbulent viscosity at cell
        symmTensor S_cell = S[celli];       // strain rate at cell

        // Compute S_ij S_ij for this cell
        scalar S_mag2 = magSqr(S_cell);

        // Effective viscosity (nu + nut)
        scalar nuEff = nu + nut_cell;

        // Dissipation formula
	scalar eps   = 2.0 * nuEff * S_mag2 ;
        epsilon_eff[celli] = eps ;

        if (eps > SMALL) 
	// Avoid division by zero
	   kolmogorov[celli] = Foam::pow( (nu*nu*nu)/eps, 0.25);
        else
	   kolmogorov[celli] = 0.0 ;	
    }

    Info << "Writing field 'epsilon_eff' and 'kolmogrov' to time " << runTime.timeName() << endl;
    epsilon_eff.write();
    kolmogorov.write();
    Info << "Done." << endl;
    return 0;
}

