    Info<< "Retrieving node positions" << endl;
    if(Pstream::master())
    {
        // --loop over nodes in the BeamDyn blade model (assumed single element)
        double posi[3], roti[3];
        for( int inode=0; inode<nnodes; ++inode ) 
        {
            // get node position
            beamDynGetNode0Position( &inode, posi, roti );

            for( int dir=0; dir<3; ++dir )
            {
                pos[inode].component(dir) = posi[dir];
                rot[inode].component(dir) = roti[dir];
            }
            r[inode] = posi[bladeDir];

            Info<< "node " << inode << " at " 
                << posi[0] << "," << posi[1] << "," << posi[2]
                << " with orientation " 
                << 180.0/pi*roti[0] 
                << "," << 180.0/pi*roti[1] 
                << "," << 180.0/pi*roti[2]
                << "  =>  r= " << r[inode]
                << endl;
        }
    }
    Pstream::scatter(pos);
    Pstream::scatter(rot);
    Pstream::scatter(r);

//
// Don't actually need any of this information...
//
//    if(Pstream::master())
//    {
//        double gllp[ngp+1], gllw[ngp+1];
//        beamDynGetGLLPts( gllp, gllw );
//        Info<< "Gauss-Lobatto-Legendre points : " << endl;
//        for( int ip=0; ip<=ngp; ++ip ) 
//        {
//            Info<< "  pt " << ip 
//                << " loc/weight : " << gllp[ip] << " " << gllw[ip] << endl;
//        }
//
//        double gp[ngp], gw[ngp];
//        beamDynGetGaussPts( gp, gw );
//        Info<< "Gauss points (for applying distributed loads) : " << endl;
//        for( int ip=0; ip<ngp; ++ip ) 
//        {
//            Info<< "  pt " << ip 
//                << " loc/weight : " << gp[ip] << " " << gw[ip] << endl;
//        }
//    }
//    Pstream::scatter(gp);
//    Pstream::scatter(gllp);

    // TODO: pass data to motion solver

