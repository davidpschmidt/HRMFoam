if (haveSeal)
{

  const scalar timeValue = U.db().time().timeOutputValue();

  Info << "Interpolated lift profile value = " << profileTable(timeValue) << endl;


// set the equilibrium value

  volScalarField eqSeal(seal) ;  // an easy construction


  if ( dimensionedScalar("triggerDisp",dimensionSet(0,0,0,0,0,0,0),profileTable(timeValue))<triggerDisp )
    {
      Info <<"Activating sealing" << endl;

      forAll(mesh.C(), celli) // for all mesh cell centres
	{

	  if( bb.containsInside(mesh.C()[celli]) )
	    {
              eqSeal.primitiveFieldRef()[celli] = 1.0;
	    }
          else
            {
              eqSeal.primitiveFieldRef()[celli] = 0.0;
            }
	}


    }
 else // lift too high-- not in the time period when squish should be activated
    {

      eqSeal = 0.0 ;

    } // end of test whether we are at such low lift, that squish should be activeted

  seal =  ( seal + eqSeal * runTime.deltaT()/slRelaxation )
           /
            (     1  + runTime.deltaT()/slRelaxation )  ;

}
