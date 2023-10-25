
Info<< "Building global boundary list" << endl;

wordList globalPatchNames;

{
#   include "Pstream.H"

    int i;
    label numPatches;
    SortableList<word> sortablePatchNames(0);
    List<wordList> gatheredPatchNames(Pstream::nProcs());

    // collect all real patches into a wordList
    wordList realPatches(0);
    forAll (mesh.boundaryMesh(), patchIndex)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchIndex];
        if (typeid(patch) != typeid(processorPolyPatch))
        {
            realPatches.setSize(realPatches.size()+1);
            realPatches[realPatches.size()-1] = patch.name();
        }
    }

    // gather real patches from each processor
    gatheredPatchNames[Pstream::myProcNo()] = realPatches;
    Pstream::gatherList(gatheredPatchNames);

    if (Pstream::master())
    {
        // store unique patches (using a HashSet)
        HashSet<word> hashSetPatchNames(0);
        forAll (gatheredPatchNames, patchListProc)
        {
            forAll (gatheredPatchNames[patchListProc], patchNameItr)
            {
                hashSetPatchNames.insert(
                    gatheredPatchNames[patchListProc][patchNameItr]
                    );
            }
        }

        // store patch names into a SortableList
        sortablePatchNames.setSize(hashSetPatchNames.size());
        i=0;
        forAllConstIters(hashSetPatchNames, iter)
        {
            sortablePatchNames[i++] = iter.key();
        }
        sortablePatchNames.sort();

        // save the number of patches
        numPatches = sortablePatchNames.size();
    }

    // scatter the number of patches to each processor
    Pstream::scatter(numPatches);

    // reset the list size on all other processors
    if (!Pstream::master())
    {
        sortablePatchNames.setSize(numPatches);
    }

    // scatter individual patch names to each processor
    // and rebuild the patch list
    word tmpPatchName;
    for (i=0; i<numPatches; i++)
    {
        if(Pstream::master())
        {
            tmpPatchName = sortablePatchNames[i];
        }
        Pstream::scatter(tmpPatchName);
        if(!Pstream::master())
        {
            sortablePatchNames[i] = tmpPatchName;
        }
    }

    // save the global patch names
    globalPatchNames = sortablePatchNames;
}

