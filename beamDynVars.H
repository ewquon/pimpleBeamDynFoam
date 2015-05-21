    vectorList pos(nnodes,vector::zero);    // current position (for calling beamDynGetNodePosition)
    vectorList rot(nnodes,vector::zero);    // current orientation (for calling beamDynGetNodePosition)
    scalarList r(nnodes,0.0);               // spanwise coordinate
    vectorList x0(nnodes,vector::zero);     // initial positions -- NOT USED?

    scalarList s(nSurfNodes, 0.0);          // parametrized position, in range [-1,1]
    double hx[nSurfNodes*nnodes];
    double htemp[nnodes];


// Unused variables:

//    scalarList gp(ngp,0.0);               // gauss point location
//    scalarList gllp(ngp,0.0);             // GLL point location (parametrized r, between [-1,1])

