    int nSurfNodes = mesh.boundaryMesh()[interfacePatchID].localPoints().size();
    //Pout<< "number of interface nodes = " << nSurfNodes << endl; // Note: will be greater than actual number of surface nodes
                                                                    // because of overlap across procs
    scalarList s(nSurfNodes, 0.0); // parametrized position, in range [-1,1]

    if( bladeR0 < 0.0 ) bladeR0 = r[0];
    if( bladeR  < 0.0 ) bladeR  = r[nnodes-1];
    scalar bladeHalfLen = (bladeR - bladeR0)/2.0;
    Info<< "Blade span : " << bladeR0 << " " << bladeR << endl;

    Info<< "NOTE: surface parametrization should be performed before the mesh is moved from its original position" << endl;
    if (nSurfNodes > 0)
    {
        const pointField& bladePoints = mesh.boundaryMesh()[interfacePatchID].localPoints();

        forAll( bladePoints, ptI )
        {
            s[ptI] = (bladePoints[ptI].component(bladeDir) - bladeR0) / bladeHalfLen - 1.0;
        }
        Pout<< "parametrized blade coordinates between [ " << min(s) << ", " << max(s) << " ] "
            << "(size=" << nSurfNodes << ")" << endl;
    }
