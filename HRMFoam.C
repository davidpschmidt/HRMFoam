/*---------------------------------------------------------------------------*\

Copyright 2006-2010 David Schmidt, Shiva Gopalakrishnan, Kshitij Neroorkar
University of Massachusetts Amherst

License
This file is part of HRMFoam, which includes software from OpenFOAM.

OpenFOAM is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OpenFOAM; if not, write to the Free Software Foundation,
Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    HRMFoam

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible flash-boiling or condensing flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "SortableList.H"
#include "./propReader/refprop/refprop.H"
#include "./propReader/utrc/utrc.H"
#include "./modelCal/modelCalc.H"
#include "turbulentFluidThermoModel.H"
#include "cellSet.H"
#include "IOobjectList.H"
#include "multivariateGaussConvectionScheme.H"
#include "processorPolyPatch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
     #include "setRootCase.H"
     #include "createTime.H"
     #include "createDynamicFvMesh.H"
     #include "readThermodynamicProperties.H"
     #include "readTransportProperties.H"
     #include "readSYconstants.H"
     #include "readSealconstants.H"
     #include "createFields.H"
     #include "setupSeal.H"

     // Read and set automatic time step
     // Changed to use variable maxCo if enabled - GLJ 7/16/18
     #include "readTimeControls_variableMaxCo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#   include "buildGlobalBoundaryList.H"
#   include "setCompressibleDeltaT.H"    // set initial timestep

#include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readHRMPISOControls.H"
        // Reread time controls from controlDict (and optionally maxCo.dat)
#       include "readTimeControls_variableMaxCo.H"
        // Calculate Courant number based only on flux divided by density
      //#       include "compressibleCourantNo.H"
#       include "fe32Courant.H"
        // Choose the largest possible dT
#       include "setDeltaT.H"

        //make phiv absolute
        fvc::makeAbsolute(phiv, U);

        // Increment time after deltaT is set
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // calculate divPhiv for topo changes with phiv from previous mesh
        divPhiv = fvc::div(phiv);

        // update mesh
	bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "seal.H"

//#       include "volContinuity.H" // AJH
#       include "continuityErrs.H" // AJH
	if(correctPhi && meshChanged)
	{
#           include"correctPhiv.H"
	}

        // make phiv relative to moving mesh
        fvc::makeRelative(phiv, U);
//        phi = phi -  fvc::interpolate(rho)*fvc::meshPhi(rho, U);

	// use a predictor/corrector step if compressible
	for ( int iStep = 0; iStep <= compressible ; iStep++)
	{

#           include "rhoEqn.H"
#           include "xEqn.H"
#           include "yEqn.H"
#           include "tracer.H"

        if(!adiabatic)
        {
#           include "hEqn.H"
        }

        //update all properties based on based on current pressure and enthalpy
        Model.update_props(h,x,y,rho);

        turbulence->correct();

        // solve for new U with lagged pressure.
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
	  + fvm::div(phi, U,"div(phi,U)")
     //     - fvm::laplacian(Model.mu(), U)
          + turbulence->divDevRhoReff(U)
          + fvm::Sp(rho*seal/drag,U)
	  // could add gravity later?
        );

        solve(UEqn == -fvc::grad(p));

	// rUA is reciporical of the diag. term of A
        rUA = 1.0/UEqn.A();

        // Pressure change for secant rule
        volScalarField dp(pSave - p);

        const dimensionedScalar pScale
        (
             "pScale",
             dimensionSet(1, -1 ,-2, 0,0,0,0),
             1.0
        );

        dimensionedScalar dpMin
        (
            //Foam::max(pScale*SMALL, (Foam::max(p) - Foam::min(p))*SMALL)
            Foam::max(pScale*Foam::sqrt(SMALL), (Foam::max(p) - Foam::min(p))*Foam::sqrt(SMALL)) // suggested by Sasan 10 Nov 2021 - see EP 1548
        );

        surfaceScalarField pFlux
        (
            fvc::interpolate(rUA)*mesh.magSf()*fvc::snGrad(p)
        );

	doubleScalar pResid = 0;
        int corr = 0;

        //make phiv absolute
        fvc::makeAbsolute(phiv, U);

	//   PISO/Secant Loop
        do
        {
	    U.storePrevIter();

            // Get U without the pressure correction
	    U = UEqn.H()*rUA;

            // Volume flux without the pressure gradient correction (absolute)
            surfaceScalarField phivStar
            (
                ( fvc::interpolate(U) & mesh.Sf())
	//	+ fvc::ddtPhiCorr(rUA,rho , U, phiv);
            );

            // update properties and calculate model term and  its derivative
            Model.update(h,rho, x, y, stabilise(dp,dpMin));
            MSave = Model.MTerm();
            dMdp =  Model.dMdp();

            // Solve pressure equation:  psi is drho/dp and omega is drho/dx

	    p.storePrevIter();

	    pFlux.storePrevIter();

	    pSave = p;

#           include "pEqn.H"

            U -= rUA*fvc::grad(p);

	    U.relax();

	    // Remove the swirl component of velocity for "wedge" cases
	    if (piso.found("removeSwirl"))
	    {
		label swirlCmpt(readLabel(piso.lookup("removeSwirl")));
		U.field().replace(swirlCmpt, 0.0);
	    }

            U.correctBoundaryConditions();

            p.relax();

	    p.correctBoundaryConditions();

            dp = pSave - p;

            #ifdef FULLDEBUG
            Info<<nl <<" Press Change: "  << Foam::max(mag(pSave - p)).value()<<endl;
            #endif
            //update all properties based on based on current pressure and enthalpy
            Model.update_props(h,x,y,rho);
            corr++;
        }
        while (pResid >= NRepsilon && corr < nCorr);

        //Make phi and phiv relative
        fvc::makeRelative(phiv,U);
        phi = fvc::interpolate(rho)*phiv;

        } // end iStep loop

        if(solveSigma)
        {
#           include "OmegaEqn.H"
        }
        Info<< nl <<"Minimum pressure is "<< Foam::min(p).value()/1.0E6
            <<" [MPa]  "
            <<      "Maximum pressure is "<< Foam::max(p).value()/1.0E6
            <<" [MPa]  "
            <<      "Minimum density is " << Foam::min(rho).value()
            <<" [kg/m3]  "
            <<      "Maximum density is " << Foam::max(rho).value()
            <<" [kg/m3]  "
            <<      "Maximum velocity is "<< Foam::max(mag(U)).value()
            <<" [m/s]";
        if(sealWithStaticMesh)
	{
		Info<< "  Maximum seal value is "<< Foam::max(seal).value();
	}
	Info << endl;

	// update mole fractions for output
	Liq_Mol_a = Model.Liq_Mol_a();

	Vap_Mol_a = Model.Vap_Mol_a();
	
        if (Foam::max(mag(U)).value() >5000.00)
        {
            Warning << "Case diverging: U > 5000 m/s - finalising" << endl;
           runTime.writeAndEnd();
        }

        runTime.write();

#       include "computeMassFlux.H"

    } // end of time stepping loop

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
