Info<< "Reading SigmaY model constants\n" << endl;

IOdictionary sigmayConstants
(
 IOobject
 (
  "sigmayConstants",
  runTime.constant(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::NO_WRITE
  )
);

dimensionedScalar Sc("Sc",dimless,VGREAT);
if (sigmayConstants.found("Sc"))
{
  Sc.read(sigmayConstants);
}

Info<< "Turbulent Schmidt number = "<< Sc.value() << endl;

Switch solveSigma("no");
solveSigma.readIfPresent("solveSigma", sigmayConstants);

dimensionedScalar Cmu("Cmu", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.09);
if (sigmayConstants.found("Cmu"))
  {
    Cmu.read(sigmayConstants);
  }

dimensionedScalar sigmaST("sigmaST", dimensionSet(1, 0, -2, 0, 0, 0, 0) ,2.89e-2);
if (sigmayConstants.found("sigmaST"))
  {
    sigmaST.read(sigmayConstants);
  }

    dimensionedScalar C4("C4", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,5.5);
if (sigmayConstants.found("C4"))
  {
    C4.read(sigmayConstants);
  }

    dimensionedScalar C5("C5", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1.0);
if (sigmayConstants.found("C5"))
  {
    C5.read(sigmayConstants);
  }

    dimensionedScalar C6("C6", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1.0);
if (sigmayConstants.found("C6"))
  {
    C6.read(sigmayConstants);
  }

// To make HRMFoam backwards compatable and avoid redeclaration of
// Sigma constants in every iteration, Sigma constants will be declared here.
// The read / write options will be conditional so that non-sigma runs don't see
// Sigma fields.

int conditionalWrite = 1;
if (! solveSigma) conditionalWrite = 0; // added DD to suppress sigma fields when solveSigma='no'

/*volScalarField Sigma
(
  IOobject
  (
    "Sigma",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0)
);*/


volScalarField Omega
(
  IOobject
  (
    "Omega",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar(dimensionSet(-1, 2, 0, 0, 0, 0, 0), Zero) // DD changed this Feb 2023
);

volScalarField Ybar
(
  IOobject
  (
    "Ybar",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), Zero)
);

volScalarField lt
(
  IOobject
  (
    "lt",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 0, 0), Zero)
);

volScalarField d32
(
  IOobject
  (
    "d32",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 0, 0), Zero),
  zeroGradientFvPatchScalarField::typeName
);

volScalarField A
(
  IOobject
  (
    "A",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), Zero)
);

volScalarField a
(
  IOobject
  (
    "a",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), Zero)
);

    /*    volScalarField rEq
    (
        IOobject
        (
            "rEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
            ),
        mesh,
        dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 0, 0), Zero)
    );*/

volScalarField sEq
(
  IOobject
  (
    "sEq",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), Zero)
);

volScalarField sMin
(
  IOobject
  (
    "sMin",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), Zero)
);

    volScalarField sInit
    (
        IOobject
        (
            "sInit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
        mesh,
        dimensionedScalar(dimensionSet(0, -1, -1, 0, 0, 0, 0), Zero)
	);

volScalarField Sigma
(
 IOobject
 (
  "Sigma",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::writeOption(conditionalWrite)
 ),
 mesh,
 dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), Zero)
 );

