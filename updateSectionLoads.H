if( Pstream::master() ) 
{

    int nnodes;
    beamDynGetNnodes(&nnodes);
    //Info<< "number of beamdyn nodes = " << nnodes << endl;

    double pos[3], rot[3];

    for( int inode=1; inode<=nnodes; ++inode ) 
    {
        //beamDynGetDispAtNode( &inode, disp, rot );
        beamDynGetNodePosition( &inode, pos, rot );
        Info<< "node " << inode << " at " 
            << pos[0] << "," << pos[1] << "," << pos[2]
            << " with orientation " 
            << 180.0/pi*rot[0] << "," << 180.0/pi*rot[1] << "," << 180.0/pi*rot[2]
            << endl;
        // TODO: specify blade axis; right now assume x
        scalar r = pos[0];
        Info<< "  using r= " << r << endl;


    }

}
