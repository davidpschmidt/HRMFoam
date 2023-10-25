{
        scalar smoothCount = 0;

        while ( (Foam::min(rho) < rhoMin) && ( smoothCount < 100) )
        {
            smoothCount++;

#ifdef FULLDEBUG
            Info << "Smoothing rho: Min = "<< Foam::min(rho).value()<<"  Max = "<< Foam::max(rho).value() << " [kg/m3] \n";
#endif
            const dimensionedScalar Dvisc("Dvisc", dimless, 1.0e-3);

            surfaceScalarField Fflag(Foam::max(scalar(1),2*fvc::interpolate(neg(rho-rhoMin))));

            surfaceScalarField rDelta2(1.0/sqr(mesh.deltaCoeffs()));

            rho += fvc::laplacian(Dvisc*rDelta2*Fflag,rho);
        }

        if (smoothCount > 0) Info << "Density smoothed in "<< smoothCount<<" iterations. \n";
        if (smoothCount == 100) Info <<"WARNING:  Density smoothing failed to sufficiently increase the minimum density. \n";
}
