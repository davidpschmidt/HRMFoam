{
    scalar netFlux = 0;
    forAll (globalPatchNames, patchI)
    {
        scalar momFlux =0;
        scalar flux = 0;
        label patchIndex =
            mesh.boundaryMesh().findPatchID(globalPatchNames[patchI]);

        if (patchIndex >= 0)
        {
            flux = sum(phi.boundaryField()[patchIndex]);
            momFlux = sum(phi.boundaryField()[patchIndex]* mag(U.boundaryField()[patchIndex]));
        }



        reduce(flux, sumOp<scalar>());
        reduce(momFlux, sumOp<scalar>());

        Info<< "Mass flux at " << globalPatchNames[patchI]
            << " = " << flux << endl;
        Info<< "Momentum flux at " << globalPatchNames[patchI]
            << " = " << momFlux << endl;

        netFlux += flux;
    }

    Info<< "Time = " << runTime.timeName()
        << " Net mass flux = " << netFlux
        << nl << endl;


    forAll (globalPatchNames, patchI)
    {
            scalar tracerMassFlux=0;
            scalar tracerValue =0;

            label patchIndex =
                    mesh.boundaryMesh().findPatchID(globalPatchNames[patchI]);

            forAll (tracerNames,tr)
            {
                    if(patchIndex >=0)
                {
                      tracerMassFlux = sum(phi.boundaryField()[patchIndex] * Tracer[tr].boundaryField()[patchIndex] );
                      tracerValue = sum(Tracer[tr].boundaryField()[patchIndex]);
                }   

                reduce(tracerMassFlux, sumOp<scalar>());
                reduce(tracerValue, sumOp<scalar>());
                
                Info<<"Mass flux through nozzle "<<Tracer[tr].name()<< " at boundary "<<globalPatchNames[patchI]<<" = "<<tracerMassFlux<<endl;                 
                Info<<"total tracer value "<<tracerValue<<endl;  
            }  

    }
  
   Info << endl;
}

