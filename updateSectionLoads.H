if( Pstream::master() ) 
{

    int nnodes;
    beamDynGetNnodes(&nnodes);
    //Info<< "number of beamdyn nodes = " << nnodes << endl;

    //dimensionedScalar pRef( dimensionedScalar(dimPressure, pRefValue) / rhoRef );
    scalar pRef( pRefValue / rhoRef );

    // setup arrays, pointers
    double pos[3], rot[3];
    double r0, r1;
    scalarList r(nnodes,0.0);
    const polyPatch& bladePatch = mesh.boundaryMesh()[interfacePatchID];
    const vectorField& bladePatchNormals = mesh.Sf().boundaryField()[interfacePatchID];

    vectorList Fp(nnodes, vector::zero);
    vectorList Fv(nnodes, vector::zero);
    vectorList Mp(nnodes, vector::zero);
    vectorList Mv(nnodes, vector::zero);

    // calculate shear stress
    Info<< "Calculating surface shear stresses" << endl;
    const volSymmTensorField Reff(turbulence->devReff());
    vectorField bladePatchShearStress = 
        (
            -mesh.Sf().boundaryField()[interfacePatchID]
            /mesh.magSf().boundaryField()[interfacePatchID]
        ) & Reff.boundaryField()[interfacePatchID];

    //
    // --loop over nodes in the BeamDyn blade model (assumed single element)
    //
    for( int inode=1; inode<=nnodes; ++inode ) 
    {
        // get node position
        beamDynGetNodePosition( &inode, pos, rot );
        Info<< "node " << inode << " at " 
            << pos[0] << "," << pos[1] << "," << pos[2]
            << " with orientation " 
            << 180.0/pi*rot[0] << "," << 180.0/pi*rot[1] << "," << 180.0/pi*rot[2]
            << endl;
        r[inode-1] = pos[bladeDir];
        Info<< "  using r= " << r << endl;
    }

    for( int inode=0; inode<nnodes; ++inode ) 
    {
        // calculate section width
        if( inode==0 )
        {
            r0 = r[0]; 
            r1 = r[1];
        }
        else if( inode==nnodes-1 )
        {
            r0 = r[nnodes-2];
            r1 = r[nnodes-1];
        }
        else
        {
            r0 = r[inode-1];
            r1 = r[inode+1];
        }
        scalar dr = (r1 - r0)/2.0;

        // 
        // --loop over faces on interface patch
        //
        forAll( bladePatch, faceI )
        {
            vector rc( bladePatch.faceCentres()[faceI] );
            if( rc[bladeDir]  >= r0 && rc[bladeDir] < r1 )
            {
                vector Sf( bladePatchNormals[faceI] ); // surface normal
                vector dm( rc - origin );

                vector dFp = rhoRef * Sf * (p.boundaryField()[interfacePatchID][faceI] - pRef) / dr;
                Fp[inode] += dFp;
                Mp[inode] += dm ^ dFp;

                vector dFv = mag(Sf) * bladePatchShearStress[faceI] / dr;
                Fv[inode] += dFv;
                Mv[inode] += dm ^ dFv;
            }

        }// end of face loop

    }

    // end loop over beamdyn nodes

}
