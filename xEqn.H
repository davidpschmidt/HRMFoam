
Model.update_props(h,x,y,rho);

{

    solve
        (
            fvm::ddt(rho, x)
          + mvConvection->fvmDiv(phi, x)

         ==

	    rho*Model.DxDt()
        );

    x = Foam::max(x,0.0);
    x = Foam::min(x,1.0);    

}

