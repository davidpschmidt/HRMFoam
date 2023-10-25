//----------------------------------------------------------------------
//  solve the Sigma (interfacial area density per unit volume) equation
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
    const volScalarField nut = turbulence->mut()/rho;
//  Turbulent length scale
    lt = Foam::pow(Cmu,3.0/4.0)*Foam::pow(turbulence->k(),3.0/2.0)/Foam::max(turbulence->epsilon(),SMALL*ulen*ulen/(utime*utime*utime));
  
//  update the Ybar and yl fields

    volScalarField yl = (1.0-y)*(1.0-x); //liquid mass fraction
    Ybar = rho*(1.0-y)*(1.0-x)/Model.rhoL(); //liquid volume fraction
    Ybar = Foam::max(Foam::min(Ybar,1.0),0.0);
    volScalarField fcntlCutOff=pos(1.0-2.004*mag(Ybar-0.5)); //cutoff for calculating the d32 
    volScalarField liqCutOff = pos(1.0-2.004*mag(yl-0.5)); //cuttoff for liquid mass fraction, to be used in the calculation of the sigma initialization term
   
    A =C4*(Cmu*fvc::grad(U)&&fvc::grad(U)*turbulence->k()/stabilise(turbulence->epsilon(),SMALL*ulen*ulen/(utime*utime*utime))); //production term due to turbulent stretching of the interface, refer [1,2]
    
    a =C5*(turbulence->epsilon()/(Foam::max(turbulence->k(),smallk))); //production term due to collisions, refer [1,2]
    
    sEq = Foam::max((C6*(Model.rhoL()+Model.rhog())*Ybar*(1.0-Ybar)*turbulence->k()/sigmaST),SMALL/ulen); //Sigma equilibrium term, refer [4]
    
    sInit = (2.0*nut*fvc::grad(yl)&&fvc::grad(yl)/Sc)*( (pos(0.001-(yl*(1.0-yl)))*6*rho*rho/(Model.rhoL()*Model.rhog()*Foam::max(lt,SMALL*ulen))) + ((1.0-pos(0.001-(yl*(1.0-yl))))*Sigma/Foam::max(yl*(1.0-yl),SMALL))); //Sigma initialization term, refer [3]

    //(a*Sigma/sEq) is the destruction term due to coalescence, refer [1,2]

//  Construct and solve the sigma transport equation
    int count = 0;
    int nLoop = 20;
    doubleScalar sigmaResid = 0;
    do{

       fvScalarMatrix sigmaEqn
        (
            fvm::ddt(Sigma)
            + fvm::div(phiv,Sigma)
            - fvm::laplacian(nut/Sc,Sigma)
	    - sInit
            - fvm::Sp((A+a),Sigma)
	    + fvm::Sp((a*Sigma/sEq),Sigma)
	    
      	 );
       sigmaResid = solve(sigmaEqn).initialResidual();
       Info << "Max Sigma: " << Foam::max(Sigma).value() << " Min Sigma: " << Foam::min(Sigma).value() << endl; 
    count++;
       }
    while(sigmaResid >= 1e-05 && count<nLoop);

    Sigma = Foam::max(Sigma,0.0/ulen); //make Sigma positive
  
//  Update equivalent D32 (SMD)
//
    d32 = (6.0*Ybar/(Foam::max(Sigma,SMALL/ulen)))*fcntlCutOff;
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
