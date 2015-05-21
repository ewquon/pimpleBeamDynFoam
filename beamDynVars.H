    vectorList pos(nnodes,vector::zero);    // current position (for calling beamDynGetNodePosition)
    vectorList rot(nnodes,vector::zero);    // current orientation (for calling beamDynGetNodePosition)
    scalarList r(nnodes,0.0);               // spanwise coordinate, used in updateSectionLoads.H

    scalarList s(BD::nSurfNodes, 0.0);          // parametrized position, in range [-1,1]
    double hx[BD::nSurfNodes*nnodes];
    double htemp[nnodes];


// Unused variables:

//    vectorList x0(nnodes,vector::zero);     // initial positions
//    scalarList gp(ngp,0.0);               // gauss point location
//    scalarList gllp(ngp,0.0);             // GLL point location (parametrized r, between [-1,1])

