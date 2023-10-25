
if (tracerON)
  {
    //for all tracers
    forAll(tracerNames,tr)
      {
	//Solve tracer equation
	fvScalarMatrix TracerEq
	  (
	   fvm::ddt(rho, Tracer[tr])
	   + mvConvection->fvmDiv(phi, Tracer[tr])
	   - fvm::laplacian(turbulence->mut()/Sc,Tracer[tr])
	   );
	solve(TracerEq,mesh.solver("Tracer"));

	Tracer[tr] = Foam::min(Tracer[tr],tracerMax);
	Tracer[tr] = Foam::max(Tracer[tr],tracerMin);


	if((cellObjects.lookup(Tracer[tr].name())) != NULL)
	  {
	    //get the current tracer cellset
	    cellSet set(*(cellObjects.lookup(Tracer[tr].name())));

	    //set the tracer value to zero for all tracers in the cellSet
	    forAll(tracerNames, trNms)
	      {

		for(cellSet::iterator cell=set.begin();cell!=set.end();cell++)
		  {
		    Tracer[trNms][cell.key()]=0.00;
		  }
	      }
	    //set the value for the correct tracer to 1
	    for(cellSet::iterator cell=set.begin();cell!=set.end();cell++)
	      {
		Tracer[tr][cell.key()]=1.00;
	      }
	  }

      }

  }
