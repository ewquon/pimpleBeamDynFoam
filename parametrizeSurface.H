    //Pout<< "number of interface nodes = " << nSurfNodes << endl; // Note: will be greater than actual number of surface nodes
                                                                    // because of overlap across procs

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

    Info<< "Retrieving element shape function" << endl;
    if (Pstream::master())
    {
        // setup initial position array -- DONT NEED THIS
//        double Nuu0[3*nnodes]; // assume 1 element
//        for( int inode=0; inode<nnodes; ++inode )
//        {
//            for( int dir=0; dir<3; ++dir )
//            {
//                Nuu0[3*inode+dir] = x0[inode][dir];
//            }
//        }

        // get shape functions
        for( int isurf=0; isurf<nSurfNodes; ++isurf )
        {
            beamDynGetShapeFunctions( &s[isurf], htemp ); 
            for( int inode=0; inode<nnodes; ++inode )
            {
                hx[isurf*nnodes + inode] = htemp[inode];
            }
        }
    }

