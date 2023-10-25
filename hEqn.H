{
  cpMix = (1-x)*(1-y)*cpL + x*(1-y)*cpV + y*Model.cpGas();

  alphaMix = Model.K()/(cpMix); //removed rho from denominator. rho IS mix rho, includes L/V/G
  //set laminar diffusivity alpha in modelCalc
  Model.set_alpha(alphaMix);


  if (U.db().foundObject<volScalarField>("alphat"))
    {
      const volScalarField& alphat = U.db().lookupObject<volScalarField>("alphat");
      alphaMix += alphat;
    }

  volTensorField gradU(fvc::grad(U));

  fvScalarMatrix hEqn
    (
     fvm::ddt(rho,h)
     + fvm::div(phi,h,"div(phi,h)")
     ==
     fvc::ddt(p)
     + fvc::div(phiv,p) - p*fvc::div(phiv+fvc::meshPhi(U))
     + fvm::laplacian(alphaMix,h) //turbulent enthalpy dispersion
     + turbulence->muEff()*((gradU + gradU.T()) && gradU) // viscous dissipation
     );

  hEqn.solve();
}

