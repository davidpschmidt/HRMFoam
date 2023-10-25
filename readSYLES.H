Info << "Reading Sigma-Y model constants" << endl;


/*Note: The values assigned to variables will treated as default. These constants can be changed in $caseDir/constant/sigmayConstants. */

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

dimensionedScalar Sc("Sc",dimless,GREAT);
if (sigmayConstants.found("Sc"))
{
  Sc = dimensionedScalar(sigmayConstants.lookup("Sc"));
}

Info << "Turbulent Schmidt number = " << Sc.value() << endl;

Switch solveSigma("no");
solveSigma.readIfPresent("solveSigma", sigmayConstants);

dimensionedScalar sigmaST("zero", dimensionSet(1, 0, -2, 0, 0, 0, 0) ,2.89e-2);
if (sigmayConstants.found("surfaceTension"))
{
  sigmaST = dimensionedScalar(sigmayConstants.lookup("surfaceTension"));
}

dimensionedScalar dnsAlpha("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,2.4);
if (sigmayConstants.found("dnsAlpha")) 
{
 dnsAlpha = dimensionedScalar(sigmayConstants.lookup("dnsAlpha"));
}

dimensionedScalar nozDia("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0), 165e-06); //Assigned the sprayG nozzle diameter as a default value
if (sigmayConstants.found("nozDia"))
{
  nozDia = dimensionedScalar(sigmayConstants.lookup("nozDia"));
}

dimensionedScalar C1("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1.0);
if (sigmayConstants.found("C1")) 
{
  C1 = dimensionedScalar(sigmayConstants.lookup("C1"));
}

dimensionedScalar C2("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1.0);
if (sigmayConstants.found("C2")) 
{
 C2 = dimensionedScalar(sigmayConstants.lookup("C2"));
}

int conditionalWrite = 1;

if (solveSigma)
{
    conditionalWrite = 0;
}

volScalarField OmegaP
(
  IOobject
  (
   "OmegaP",
   runTime.timeName(),
   mesh,
   IOobject::READ_IF_PRESENT,
   IOobject::writeOption(conditionalWrite)
  ),
  mesh,
  dimensionedScalar("zero", dimensionSet(-1, 2, 0, 0, 0, 0, 0) ,0.0)
);

volScalarField Ybar
(
 IOobject
 (
  "Ybar",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::writeOption(conditionalWrite)
 ),
 mesh,
 dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0)
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
 dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0) ,0.0),
 zeroGradientFvPatchScalarField::typeName
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
 dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
);

volScalarField SigmaEq
(
 IOobject
 (
  "SigmaEq",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0) ,0.0)
);

volScalarField CSigma
(
 IOobject
 (
  "CSigma",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
);

volScalarField OmegaMin
(
 IOobject
 (
  "OmegaMin",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::writeOption(conditionalWrite)
 ),
 mesh,
 dimensionedScalar("zero", dimensionSet(-1, 2, 0, 0, 0, 0, 0) ,0.0),
 calculatedFvPatchScalarField::typeName
);

volScalarField nD
(
 IOobject
 (
  "nD",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::writeOption(conditionalWrite)
 ),
 mesh,
 dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0),
 zeroGradientFvPatchScalarField::typeName
);

volScalarField LESdelta
(
 IOobject
 (
  "LESdelta",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("LESdelta", dimLength, SMALL),
 zeroGradientFvPatchScalarField::typeName
);

volScalarField cellVolume
(
 IOobject
 (
  "cellVolume",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("small",dimVolume,SMALL),
 calculatedFvPatchScalarField::typeName
);

volScalarField cutOff
(
 IOobject
 (
  "cutOff",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
 zeroGradientFvPatchScalarField::typeName
);
