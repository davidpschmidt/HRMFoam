    const dictionary piso = mesh.solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr = 0;
    if (piso.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(piso.lookup("nNonOrthogonalCorrectors"));
    }

    doubleScalar NRepsilon = 1.0e-6;
    if (piso.found("NRepsilon"))
    {
        NRepsilon = (piso.get<doubleScalar>("NRepsilon"));
    }

    if (piso.found("removeSwirl"))
    {
	label swirlCmpt(readLabel(piso.lookup("removeSwirl")));
	Info<< "Removing swirl component-" << swirlCmpt << " of U" << endl;

    }

    bool correctPhi = false;
    if (piso.found("correctPhi"))
      {
         correctPhi = Switch(piso.lookup("correctPhi"));
      }
