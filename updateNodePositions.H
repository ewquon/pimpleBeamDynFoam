    Info<< "Retrieving node positions" << endl;
    if(Pstream::master())
    {
        // --loop over nodes in the BeamDyn blade model (assumed single element)
        double pos[3], rot[3];
        for( int inode=1; inode<=nnodes; ++inode ) 
        {
            // get node position
            beamDynGetNodePosition( &inode, pos, rot );
            Info<< "node " << inode-1 << " at " 
                << pos[0] << "," << pos[1] << "," << pos[2]
                << " with orientation " 
                << 180.0/pi*rot[0] << "," << 180.0/pi*rot[1] << "," << 180.0/pi*rot[2]
                << endl;
            r[inode-1] = pos[bladeDir];
            Info<< "  using r= " << r[inode-1] << endl;
        }
    }
    Pstream::scatter(r);

    if(Pstream::master())
    {
        double gllp[ngp+1], gllw[ngp+1];
        beamDynGetGLLPts( gllp, gllw );
        Info<< "Gauss-Lobatto-Legendre points : " << endl;
        for( int ip=0; ip<=ngp; ++ip ) 
        {
            Info<< "  pt " << ip 
                << " loc/weight : " << gllp[ip] << " " << gllw[ip] << endl;
        }

        double gp[ngp], gw[ngp];
        beamDynGetGaussPts( gp, gw );
        Info<< "Gauss points (for applying distributed loads) : " << endl;
        for( int ip=0; ip<ngp; ++ip ) 
        {
            Info<< "  pt " << ip 
                << " loc/weight : " << gp[ip] << " " << gw[ip] << endl;
        }
    }
    Pstream::scatter(gp);
    Pstream::scatter(gllp);

