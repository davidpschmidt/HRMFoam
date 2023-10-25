ScLimited = Foam::max( turbulence->mut()*Foam::max(mag(fvc::grad(y)),gyMin) /
                                Foam::max(
                                           sqrt(2*turbulence->k()/3) * y*(1-y) , uMin)
                                         /rho
                                , Sc);

// Iteratively correct for non-orthogonality
for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
{
  // this is the partial of rho with respect to y divided by rho squared
  volScalarField DrhoDyOverRhoSq(-(1/Model.rhog() - x/Model.rhoV() + (x-1)/Model.rhoL() ));

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p)
    );

    if(compressible)
    {
     rhoEOS = 1/(  (x/Model.rhoV() + (1-x)/Model.rhoL())*(1-y) + y/Model.rhog() );

       // update psi in modelCalc
      //Model.update_psi(rhoEOS,y,x,h);
      Model.update_psi(rho,y,x,h); // suggested by Sasan 21 Oct 2021 - see EP 1548

      //uncomment psi2phase here and in createFields to have psi written for post-processing
      //psi2phase=Model.psi();

      surfaceScalarField phivRel = phiv;
      fvc::makeRelative(phivRel,U);

    // This is the pressure equation with compressibility
	pResid = solve
        (
	 (Model.psi()/rho)*
           (
	     (fvm::ddt(p) + fvm::div(phivRel, p) - fvm::SuSp(fvc::div(phiv), p))
           )

	 + fvc::div(phivStar)
	 - pEqn
	 + MSave
	 + fvm::SuSp(dMdp, p)
         - dMdp*p
	 + DrhoDyOverRhoSq*fvc::laplacian(turbulence->mut()/ScLimited,y)*Foam::pos(DrhoDyOverRhoSq*fvc::laplacian(turbulence->mut()/ScLimited,y)*(rhoEOS-rho))
        ).initialResidual();
    }
    else
    {
	// Incompressible pressure equation.
	pResid = solve
	(
	   fvc::div(phivStar)
	 - pEqn
	 + MSave
	 + fvm::Sp(dMdp,p)
	 - dMdp*p
	 + DrhoDyOverRhoSq*fvc::laplacian(turbulence->mut()/ScLimited,y)
	).initialResidual();
    }

    // On the last non-orthogonality correction, correct
    // the flux using the most up-to-date pressure
    if (nonOrth == nNonOrthCorr)
    {
	pFlux = pEqn.flux();
	pFlux.relax();
	phiv = phivStar - pFlux;
    }
}

