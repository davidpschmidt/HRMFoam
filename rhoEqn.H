
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phiv,
        mesh.divScheme("div(phiv,rho)")
    )
);


        // Solve continuity for updated density

        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho)
          + mvConvection->fvmDiv(phiv,rho)
        );


        solve(rhoEqn);


        phi = rhoEqn.flux();


        if(lBoundRho)
        {
#          include "lBoundRho.H"
        }

        // This will ensure that density remains positive.
        // This does not conserve mass
        rho == Foam::max(rho, rhoMin);

        //Some bcs use internal field values for calculating boundary values
        //these need to be updated

        rho.correctBoundaryConditions();  //will it be a problem that phi's are relative? -DPS
