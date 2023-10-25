//----------------------------------------------------------------------
//  solve the Omega (interfacial area density per unit volume) equation
//----------------------------------------------------------------------

{

//  trivially small constant used as lower bound
    const dimensionedScalar smallk
        (
            "smallk",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            scalar(SMALL)
            );

//
//  some unit quantities
//
    const dimensionedScalar ulen  //For divide by zero
        (
            "ulen",
            dimensionSet(0,1,0,0,0,0,0),
            scalar(1.0)
            );
    const dimensionedScalar utime  //For divide by zero
        (
            "utime",
            dimensionSet(0,0,1,0,0,0,0),
            scalar(1.0)
            );

//  Turbulent kinematic viscosity
    const volScalarField nut(turbulence->mut()/rho);
//  Turbulent length scale
    lt = Foam::pow(Cmu,3.0/4.0)*Foam::pow(turbulence->k(),3.0/2.0)/Foam::max(turbulence->epsilon(),SMALL*ulen*ulen/(utime*utime*utime));

//  update the Ybar and yl fields

    volScalarField yl((1.0-y)*(1.0-x)); //liquid mass fraction
    Ybar = rho*(1.0-y)*(1.0-x)/Model.rhoL(); //liquid volume fraction
    Ybar = Foam::max(Foam::min(Ybar,1.0),0.0);

    // these cutoff functions = 1 when the Y value is from 0.001 to 0.999
    volScalarField fcntlCutOff(pos(1.0-2.004*mag(Ybar-0.5))); //cutoff for calculating Sigma
    volScalarField liqCutOff(pos(1.0-2.004*mag(yl-0.5))); //cuttoff for liquid mass fraction, to be used in the calculation of the Omega initialization term

    //production term due to turbulent stretching of the interface, refer [1,2]
    A =C4*(Cmu*fvc::grad(U)&&fvc::grad(U)*turbulence->k()/stabilise(turbulence->epsilon(),SMALL*ulen*ulen/(utime*utime*utime)));

    //production term due to collisions, refer [1,2]
    a =C5*(turbulence->epsilon()/(Foam::max(turbulence->k(),smallk)));

    //Omega equilibrium term, refer [4]
    sEq = Foam::max((C6*(Model.rhoL()+Model.rhog())*Ybar*(1.0-Ybar)*turbulence->k()/sigmaST),SMALL/ulen);

    //Omega initialization term, refer [3]
    sInit = (2.0*rho*nut*fvc::grad(yl)&&fvc::grad(yl)/Sc)*
      (
              pos(0.001-(yl*(1.0-yl)))      *    6*rho/(Model.rhoL()*Model.rhog()*Foam::max(lt,SMALL*ulen))
       + (1.0-pos(0.001-(yl*(1.0-yl))))     *    Omega/Foam::max(yl*(1.0-yl),SMALL)
      );

    //(a*Omega/sEq) is the destruction term due to coalescence, refer [1,2]

//  Construct and solve the Omega transport equation
    int count = 0;
    int nLoop = 4;
    doubleScalar OmegaResid = 0;

    const dimensionedScalar OmegaMin  //For divide by zero
        (
            "OmegaMin",
            dimensionSet(-1, 2, 0 , 0, 0, 0, 0),
            scalar(1)
            );

    // lets enter the equation solution process with a positive Omega so that the non-linear term has a correct sign
    Omega = Foam::max(Omega,OmegaMin);

    Omega.storePrevIter();

    do{

        fvScalarMatrix OmegaEqn
        (
	    fvm::ddt(rho,Omega)
            + fvm::div(phi,Omega)
            - fvm::laplacian(turbulence->mut()/Sc,Omega)
	    - sInit
            - fvm::Sp((A+a)*rho,Omega)
	    + fvm::Sp((a*Omega*sqr(rho)/sEq),Omega)

      	 );

       OmegaResid = solve(OmegaEqn).initialResidual();

       Omega.relax();

       count++;
    }
    while(OmegaResid >= 1e-02 && count<nLoop);

    Info << "Max Omega: " << Foam::max(Omega).value() << " Min Omega: " << Foam::min(Omega).value() << endl;

    //correction to prevent upstream prediction of Omega or Sigma, and also make Omega positive
    Sigma = (Foam::max(Omega,OmegaMin)*rho)*fcntlCutOff;

//  Update equivalent D32 (SMD)
//

    d32 = (6.0*Ybar/(rho*Foam::max(Omega,OmegaMin)))*fcntlCutOff;

    forAll(mesh.C(),celli)
       {
	 if (Ybar[celli]<1e-04)
        {
            d32[celli]=0;
        }
	}

    // References:
    // [1] Vallet A, Burluka A A and Borghi R, "Development of a Eulerian model for the "atomization" of a liquid jet", AAS, 11, pp. 619-642, 2001.
    // [2] Demoulin, Beau, et al., "A new model for turbulent flows with large density fluctuations: application to liquid atomization", AAS, 17, pp. 315-345, 2007.
    // [3] Beau, Menard, Lebas, et al., "Numerical jet atomization. Part II: Modeling information and comparison with DNS results", joint US-Eurpoean Fluids Engineering Summer Meeting (FEDSM), FEDSM2006-98166, ASME Proceedings, July 17-20, Miami, FL, 2006.
    // [4] Duret B, Revillon J, Menard T, Demoulin F X, "Improving primary atomization modeling through DNS of two-phase flows", Int. J. Multiphase Flow, 55, pp. 130-137, 2013.
}
