    vectorList pos(BD::nnodes,vector::zero);    // current position (for calling beamDynGetNodePosition)
    vectorList rot(BD::nnodes,vector::zero);    // current orientation (for calling beamDynGetNodePosition)
    scalarList r(BD::nnodes,0.0);               // spanwise coordinate, used in updateSectionLoads.H

    BD::nSurfNodes = mesh.boundaryMesh()[BD::interfacePatchID].localPoints().size();

    scalarList s(BD::nSurfNodes, 0.0);          // parametrized position, in range [-1,1]
    double hx[BD::nSurfNodes*BD::nnodes];       // shape function array
    double htemp[BD::nnodes];                   // temporary shape function variable


// Unused variables:

//    vectorList x0(nnodes,vector::zero);     // initial positions
//    scalarList gp(ngp,0.0);               // gauss point location
//    scalarList gllp(ngp,0.0);             // GLL point location (parametrized r, between [-1,1])

