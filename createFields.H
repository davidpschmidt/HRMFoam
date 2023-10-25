    //Make conditional write option so that HRMFoam doesn't write so much data to file
    Switch writeFewerFields("no");
    writeFewerFields.readIfPresent("writeFewerFields", runTime.controlDict());

    conditionalWrite = IOobject::AUTO_WRITE;

    if (writeFewerFields)
    {
      conditionalWrite = IOobject::NO_WRITE;
    }

    // x is read from an input file
    Info<< "Reading field x\n" << endl;
    volScalarField x
    (
        IOobject
        (
            "x",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // y is read from an input file
    Info<< "Reading field y\n" << endl;
    volScalarField y
    (
        IOobject
        (
            "y",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // velocity is read from an input file
    Info<< "Reading field U\n" << endl;
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

    // enthalpy
    Info<< "Reading field h\n" << endl;
    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

const dimensionedScalar cp_fake ("cp_fake", dimensionSet(0, 2, -2, -1, 0, 0, 0) ,1.0E3);

{
Info<< "Writing initial dummy T field\n" << endl;
 volScalarField T
   (
      IOobject
      (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
       ),
       h/cp_fake

     );

    T.write();
}
    // density
    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Model term
    volScalarField MSave
    (
        IOobject
        (
            "MSave",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0) ,0.0)
    );

    //derivative of the model term
    volScalarField dMdp
    (
        IOobject
        (
            "dMdp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::writeOption(conditionalWrite)
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 1, 1, 0, 0, 0, 0) ,0.0)
    );


    // volume flux with pressure correction
#   include "createPhiv.H"
    // create the initial flux field phi with variable density taken into account
    // creates surface scalar field by reading if present, or calculating from rho
    // and U using interpolation
#   include "compressibleCreatePhi.H"

    //object of class modelCalc to calculate the model term
    modelCalc Model
    (
        theta_0,
        pNonDimExp,
        pcritical,
        compressible,
        theta_floor,
        alphaFracExp,
        h,
        rho,
        x,
	y,
        U,
	conditionalWrite
    );

Info << "Creating pressure field \n";

    volScalarField& p = Model.p();

    // used for calculating change in pressure
    volScalarField pSave("pSave",p);

//update all properties based on current pressure and enthalpy
//has to be here because the turbulence model needs the vicosity to be updated
//the pressure is not needed to be sent in since modelCalc has access to the current pressure
Info << "Updating properties \n";
Model.update_props(h,x,y,rho);
Model.update_psi(rho,y,x,h);
//set alpha in modelCalc
volScalarField cpMix((1-x)*(1-y)*cpL + x*(1-y)*cpV + y*Model.cpGas());
volScalarField alphaMix(Model.K()/(cpMix)); //removed rho from denominator. rho IS mix rho, includes L/V/G
Model.set_alpha(alphaMix);

/*volScalarField psi2phase
(
    IOobject
    (
        "psi2phase",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Model.psi()
);*/

    Info<< "Creating turbulence model\n" << endl;

    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            Model
        )
    );

   Info<<"done creating turbulence model "<<endl;




    volScalarField rUA
    (
      IOobject
      (
        "rUA",
        runTime.timeName(),
        mesh
       ),
      mesh,
      dimensionedScalar("rUA",dimensionSet(-1,3,1,0,0,0,0),1.0),
      zeroGradientFvPatchScalarField::typeName
     );


 volScalarField divPhiv
    (
      IOobject
      (
        "divPhiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
       ),
      fvc::div(phiv)
     );


//this is for generating a list of cellSets for the tracers
    IOobjectList objects
    (
           mesh,
           mesh.facesInstance(),
           polyMesh::meshSubDir/"sets"
     );

    IOobjectList cellObjects(objects.lookupClass(cellSet::typeName));




   PtrList<volScalarField>Tracer(tracerNames.size());
   const dimensionedScalar tracerMin("tracerMin", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0);
   const dimensionedScalar tracerMax("tracerMax", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,1);


//defining a list of tracerss
   forAll (tracerNames, tr)
   {

       IOobject trheader
       (
            tracerNames[tr],
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );
//if the names of tracers  can be found in the time directory (i.e. 0 for start)
//then read from there

        if (trheader.typeHeaderOk<volScalarField>(true))
        {
            Tracer.set(tr, new volScalarField(trheader, mesh));
        }
//else create the new variables and initialize them to 1 given in
        else
        {
            Tracer.set
            (
                tr,
                new volScalarField
                (
                    IOobject
                    (
                        tracerNames[tr],
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "s",
                        dimensionSet(0,0,0,0,0,0,0),
                        1.00
                    ),
//this has to be a word list of boundary names
                    zeroGradientFvPatchVectorField::typeName
                )
            );

            // Use initOrientations_[oI] to initialise Q_[oI]
        }

   }


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
fields.add(rho);
fields.add(x);
fields.add(y);
   forAll (tracerNames, tr)
   {
     fields.add(Tracer[tr]);
   }



const dimensionedScalar gyMin
(
"gyMin",
dimensionSet(0,-1,0,0,0,0,0),
SMALL
);

const dimensionedScalar uMin
(
"uMin",
dimensionSet(0,1,-1,0,0,0,0),
SMALL
);
volScalarField ScLimited
   (
      IOobject
      (
        "ScLimited",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::writeOption(conditionalWrite)
       ),
       Foam::max( turbulence->mut()*mag(fvc::grad(y)) /
                               ( Foam::max( sqrt(2*turbulence->k()/3) * y*(1-y) ), uMin) /rho
                                , Sc)
     );

    volScalarField rhoEOS
    (
        IOobject
        (
            "rhoEOS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
1/(  (x/Model.rhoV() + (1-x)/Model.rhoL())*(1-y) + y/Model.rhog() )
    );
p.correctBoundaryConditions(); //Fix crash on first time step for high core counts

// Force reading of g even if not used
const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);
Info << "Gravity: " << g << " (should be zero)" << endl;

// mole fraction of species a in the liquid phase

volScalarField Liq_Mol_a
    (
      IOobject
      (
        "Liq_Mol_a",
		runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::writeOption(conditionalWrite)
       ),
      mesh,
      dimensionedScalar("Liq_Mol_a",dimensionSet(0,0,0,0,0,0,0),0.0),
      zeroGradientFvPatchScalarField::typeName
     );



// mole fraction of species a in the liquid phase

volScalarField Vap_Mol_a
    (
      IOobject
      (
        "Vap_Mol_a",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::writeOption(conditionalWrite)
       ),
      mesh,
      dimensionedScalar("Vap_Mol_a",dimensionSet(0,0,0,0,0,0,0),0.0),
      zeroGradientFvPatchScalarField::typeName
     );

