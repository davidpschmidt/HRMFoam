Info<< "Reading thermophysicalProperties\n" << endl;

IOdictionary thermophysicalProperties
(
 IOobject
 (
  "thermophysicalProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );


dimensionedScalar theta_0("theta_0", thermophysicalProperties);

//variables required for the D-Zcorrelation
dimensionedScalar alphaFracExp("alphaFracExp", thermophysicalProperties);

dimensionedScalar pNonDimExp("pNonDimExp", thermophysicalProperties);

//if pcritical is yes, then the value will be read in from the fluid.dat
Switch pcritical
(
 thermophysicalProperties.lookup("pcritical")
 );


Switch tracerON("no");
tracerON.readIfPresent("tracerON", thermophysicalProperties) ;

wordList tracerNames(thermophysicalProperties.lookup("tracerNames"));


//if adiabatic is yes, the enthalpy equation if not solved
Switch adiabatic
(
 thermophysicalProperties.lookup("adiabatic")
 );

//if lBoundRho is yes, the density will be corrected using the mass conserving
//procedure, see lBoundRho.H
//otherwise just max(rho,rhoMin) will be used
Switch lBoundRho("lBoundRho", thermophysicalProperties);

//done with variables required for D-Z correlation
dimensionedScalar theta_floor("theta_floor", thermophysicalProperties);

dimensionedScalar rhoMin("rhoMin", thermophysicalProperties);

Switch compressible
(
 thermophysicalProperties.lookup("compressible")
 );

dimensionedScalar cpL("cpL",dimSpecificHeatCapacity,VGREAT);

dimensionedScalar cpV("cpV",dimSpecificHeatCapacity,VGREAT);

if (!adiabatic)
  {
    //hardcode cpL cpV, for now -> move these to reading from fluid.dat when not garbage
    Info <<"Enthalpy equation activated (not adiabatic).  Requiring cpL and cpV inputs. \n";

    cpL.read(thermophysicalProperties);

    cpV.read(thermophysicalProperties);

  }

dimensionedScalar cpGas("cpGas", thermophysicalProperties);


