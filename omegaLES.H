/* Calulation of the mass based liquid-gas interface-area-density (Omega) in LES framework */

{

// declaration of small/unit quantities needed to avoid segmentation fault(s)
const dimensionedScalar smallk
  (
   "smallk",
   dimensionSet(0, 2, -2, 0, 0, 0, 0),
   scalar(SMALL)
  );

const dimensionedScalar ulen
  (
   "ulen",
   dimensionSet(0, 1, 0, 0, 0, 0, 0),
   scalar(1.0)
  );

const dimensionedScalar utime
  (
   "utime",
   dimensionSet(0, 0, 1, 0, 0, 0, 0),
   scalar(1.0)
  );
const dimensionedScalar smallYbar
  (
   "smallYbar",
   dimensionSet(0, 0, 0, 0, 0, 0, 0),
   scalar(SMALL)
  );

 //Info << "OmegaP Start" << endl;

// calulating the LES delta value

  if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
  {
      const volScalarField& delta = mesh.objectRegistry::lookupObject<volScalarField>("delta");
      LESdelta.internalField() = delta.internalField();
      LESdelta.correctBoundaryConditions();
  }  

  //Info << "LESdelta" << endl;

  //Info << "Max. LESdelta " << Foam::max(LESdelta).value() << " Min. LESdelta " << Foam::min(LESdelta).value() << endl;

// calulating the liquid volume fraction (Ybar)
  Ybar = rhoEOS*(1.0-y)*(1.0-x)/Model.rhoL();  
   
// bounding Ybar between [0.0,1.0]
  Ybar = Foam::max(Foam::min(Ybar,1.0),0.0);
  Ybar.storeOldTime();

volScalarField yl = (1.0-y)*(1.0-x); //liquid mass fraction

  //Info << "Ybar" << endl;

// calculating the minimum value of Omega 
  OmegaMin = ((dnsAlpha/rho)*sqrt(Ybar*(1.0-Ybar)))/LESdelta; //Foam::max(LESdelta,SMALL*ulen);
  OmegaMin.correctBoundaryConditions();

  Info << "Max. OmegaMin " << Foam::max(OmegaMin).value() << " Min. OmegaMin " << Foam::min(OmegaMin).value() << endl;

// calculation of the turbulent time scale using CSigma
  CSigma = C1 * turbulence->epsilon()/Foam::max(turbulence->k(),smallk); 

  //Info << "CSigma" << endl;

//calculation of the Sigma Equilibrium
  SigmaEq = OmegaMin*rho + Foam::max((C2*(Model.rhoL()+Model.rhog())*Ybar*(1.0-Ybar)*turbulence->k()/sigmaST),SMALL/ulen);

//Info << "SigmaEq" << endl;

//OmegaP equation

int count = 0;
int nLoop = 5;
double OmegaPResid = 0.0;
OmegaP.storePrevIter();

do {
  fvScalarMatrix OmegaPEqn
  (

    fvm::ddt(rho,OmegaP)
    + fvm::div(phi,OmegaP)
    - fvm::laplacian(turbulence->mut()/ScLimited,OmegaP)
    + fvm::Sp(2*CSigma*rho*rho*OmegaMin/SigmaEq,OmegaP)
    - fvm::Sp(CSigma*rho,OmegaP)
    + fvm::Sp(CSigma*rho*rho*OmegaP/SigmaEq,OmegaP)
    - CSigma*OmegaMin*rho*(1.0-(OmegaMin*rho/SigmaEq))
  );

  //solve(OmegaPEqn==-(Ys&fvc::grad((rho/Foam::max(Ybar,smallYbar))*OmegaP)));
  //solve(OmegaPEqn);

  //Info << "OmegaPEqn" << endl;
  OmegaPResid = solve(OmegaPEqn).initialResidual();

  OmegaP.relax();

  count++;
}
while(OmegaPResid >= 1e-03 && count<nLoop);

Info << "Maximum OmegaP: " << Foam::max(OmegaP).value() << " Minimum OmegaP: " << Foam::min(OmegaP).value() << endl; 
   
//unit quantities of Omega and Sigma for avoiding division by zero

const dimensionedScalar unitOmega
  (
   "unitOmega",
   dimensionSet(-1, 2, 0, 0, 0, 0, 0),
   scalar(1.0)
  );

const dimensionedScalar unitSigma
  (
   "unitSigma",
   dimensionSet(0, -1, 0, 0, 0, 0, 0),
   scalar(1.0)
  );

  volScalarField cutOff = pos(1.0 - 2.005*mag(Ybar-0.5));
  volScalarField liqCutOff = pos(1.0 - 2.005*mag(yl-0.5));
  Sigma = Foam::max(((OmegaMin + OmegaP)*rho),OmegaMin*rho) * cutOff; //making the total Sigma positive and zero in the case of pure liquid or pure gas

  Info << "Max. Sigma " << Foam::max(Sigma).value() << " Min. Sigma " << Foam::min(Sigma).value() << endl;

// calculation of D32
  d32 = 6*Ybar/(Foam::max(Sigma,unitSigma));
  d32 = Foam::min(d32,1.89*nozDia) * cutOff; //limiting maximum value of d32 to the Raleigh diam., and to value of 0 in the case pure liquid or gas.
  //Info << "d32" << endl;

// calculating the number of droplets (nD)
  cellVolume.internalField() = mesh.V();
  cellVolume.correctBoundaryConditions();
  nD = (1.0/mathematicalConstant::pi)*(6.0*Ybar*cellVolume/Foam::max(d32*d32*d32,SMALL*ulen*ulen*ulen)) * cutOff; //zero droplets in the case of pure liquid or gas
  //Info << "nD" << endl;

}
