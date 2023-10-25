
ScLimited = Foam::max( turbulence->mut()*Foam::max(mag(fvc::grad(y)),gyMin) /
                                Foam::max(
                                           sqrt(2*turbulence->k()/3) * y*(1-y) , uMin
                                         )
                                 /rho
                                , Sc);
{

  volScalarField Dcoeff(turbulence->mut()/ScLimited);

fvScalarMatrix yEqn
(
 fvm::ddt(rho,y)
 + mvConvection->fvmDiv(phi, y)
 - fvm::laplacian(Dcoeff,y)
 );


solve(yEqn);


y=Foam::max(y,0.0);

y=Foam::min(y,1.0);
}
