/*---------------------------------------------------------------------------*\

Global
    setCompressibleDeltaT

Description
    Set the initial timestep corresponding to the timestep adjustment
    algorithm in setCompressibleDeltaT

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
#   include "compressibleCourantNo.H"

    if (CoNum > SMALL)
    {
        runTime.setDeltaT
        (
            Foam::min
            (
                maxCo*runTime.deltaT().value()/CoNum,
                maxDeltaT
            )
        );
    }
}

// ************************************************************************* //
