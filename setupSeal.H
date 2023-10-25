volScalarField seal
 (
     IOobject
     (
         "seal",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0) ,0.0),
     zeroGradientFvPatchVectorField::typeName
 );

const boundBox bb(lowerCorner,upperCorner);//creating the sealing box

IOdictionary dMDict
(
IOobject
(
"dynamicMeshDict",
runTime.constant(),
mesh,
IOobject::MUST_READ,
IOobject::NO_WRITE
)
);

//Only set up sealing if seal constants file found
interpolationTable<scalar> profileTable;
if (haveSeal){
word motionType = dMDict.get<word>("dynamicFvMesh");
if(motionType == word("dynamicMotionSolverFvMesh")){
 Info << "Laplace detected.  Setting up sealing. " << endl;
 Info << "Converting NeedleVelocity.dat to displacement." << endl;
 Info << "To change table name or axis of lift, see setupSeal.H" << endl;

 fileName profileName("NeedleVelocity.dat");
 interpolationTable<vector> velocityTable(profileName);
 //Loop through the velocity table, keeping a running total of displacement
 //for the lift axis.  Write this displacement to a new scalar interpolationTable
 //for each time in the original table.
 List<Tuple2<scalar, scalar> > workingList(label(1));
 List<Tuple2<scalar,scalar> > theList;
 //First one is 0,0 to be safe
 Tuple2<scalar,scalar> theTuple(scalar(0),scalar(0));
 workingList[0] = theTuple;
 theList.append(workingList);
 scalar theTime(0);//initialization??--cm
 scalar lastTime(0);//initialization??--cm
 scalar dispTrack(0);//initialization??--cm
 for (label i=1; i<velocityTable.size(); i++){
     if (i!=1) {lastTime = theTime;}
     theTime = velocityTable[i].first();//Extract the time data from the velocity profile table
     dispTrack = dispTrack + (velocityTable[i].second().y()*(theTime-lastTime));//Extract the Y-Axis Velocity
     //Uncomment to print out calculated displacements
     //Info << dispTrack <<endl;
     workingList[0]=Tuple2<scalar,scalar> (theTime,dispTrack);//List created with time and displacement
     theList.append(workingList);//Saved to the main list
 }
 profileTable = interpolationTable<scalar>(theList,bounds::repeatableBounding::WARN,fileName("Derp"));//Derp is just a dummy name
 profileTable.check();
}

else if (motionType == word("pacManTopoFvMesh3")){
 Info << "PacMan detected.  Setting up sealing." << endl;
 dictionary pacManDict
 (
         dMDict.subOrEmptyDict("pacManTopoFvMesh3Coeffs")
 );
 dictionary needleSubDict
 (
         pacManDict.subOrEmptyDict("needleMotionCoeffs")
 );


 fileName profileName;
 profileName = static_cast<Foam::fileName>(needleSubDict.lookup("profileFile"));
 Info << "Opening table file " << profileName << endl;
 interpolationTable<scalar> localTable(profileName);
 profileTable = localTable;
}

// DD added support for 'sealWithStaticMesh' to use seal as drag for initial conditions.
else if ((motionType == word("staticFvMesh")) && sealWithStaticMesh){
 Info << "Static mesh with sealing applied." << endl;
 // Start with sealing active for the first few timesteps, then permit relaxation using slRelaxation
 List<Tuple2<scalar,scalar> > theList;
 Tuple2<scalar,scalar> activateSeal(scalar(0),scalar(-1));
 theList.append(activateSeal);
 Tuple2<scalar,scalar> releaseSeal(runTime.deltaT().value()*10, scalar(0));
 theList.append(releaseSeal);
 profileTable = interpolationTable<scalar>(theList,bounds::repeatableBounding::CLAMP,fileName("Foo"));
}

else {
 Info << "Invalid motion detected.  Sealing not applied." << endl;
 Info << motionType << endl;
 Info << sealWithStaticMesh << endl;
 haveSeal = false;
}


// initialize seal
if (haveSeal){
const scalar timeValue = U.db().time().timeOutputValue();

Info << "Interpolated lift profile value = " << profileTable(timeValue) << endl;

if ( dimensionedScalar("triggerDisp",dimensionSet(0,0,0,0,0,0,0),profileTable(timeValue)) < triggerDisp )
  {
    Info <<"Activating seal" << endl;

    forAll(mesh.C(), celli) // for all mesh cell centres
      {

  scalar targetValue = 0;

  if( bb.containsInside(mesh.C()[celli]) )
    {
      targetValue = 1;
    }

  seal.primitiveFieldRef()[celli] = targetValue ;

      }

  }
else // lift too high-- not in the time period when seal should be activated
  {

    forAll(mesh.C(), celli) // for all mesh cell centres
      {

  scalar targetValue = 0;  // If above sealing lift, target always zero

  seal.primitiveFieldRef()[celli] = targetValue ;

      }

  } // end of test whether we are at such low lift, that seal should be activeted
}
// note faces that border a seal cell, possibly lying between a seal and non-seal cell
} //end of first haveSeal check
