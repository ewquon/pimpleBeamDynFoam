    scalar pRef( pRefValue / rhoRef );

    // setup arrays, pointers
    double r0, r1;
    const polyPatch& bladePatch = mesh.boundaryMesh()[interfacePatchID];
    const vectorField& bladePatchNormals = mesh.Sf().boundaryField()[interfacePatchID];

    // calculate shear stress
    //   note: devReff returns the effective stress tensor including the laminar stress
    //   note: face normals point _outside_ the computational domain
    Info<< "Calculating surface shear stresses" << endl;
    const volSymmTensorField Reff(turbulence->devReff());
//    vectorField bladePatchShearStress = 
//        (
//            -mesh.Sf().boundaryField()[interfacePatchID]
//            /mesh.magSf().boundaryField()[interfacePatchID]
//        ) & Reff.boundaryField()[interfacePatchID];

// normal points into surface, i.e. the direction the fluid is pushing on the wall
// this matches the 'forces' function object implementation
// also, no need to normalize by magSf since we multiply by mag(Sf) later
    vectorField bladePatchShearStress = 
        mesh.Sf().boundaryField()[interfacePatchID]
        & Reff.boundaryField()[interfacePatchID];

    //
    // --loop over nodes in the BeamDyn blade model, assumed single element
    //   i.e., nnodes = nodes_elem = order_elem+1 = ngp+1
    //
    Info<< "Integrating sectional loads" << endl;
    for( int ig=0; ig<nnodes-1; ++ig ) 
    {
        vector Fp(vector::zero);
        vector Fv(vector::zero);
        vector Mp(vector::zero);
        vector Mv(vector::zero);

        r0 = r[ig];
        r1 = r[ig+1];
        scalar dr = r1 - r0; // note: this is the width of the integration segment 
                             //       corresponding to the Gauss pt

        // 
        // --loop over faces on interface patch
        //
        int nFacesFound = 0;
        forAll( bladePatch, faceI )
        {
            vector rc( bladePatch.faceCentres()[faceI] );
            if( rc[bladeDir] >= r0 && rc[bladeDir] < r1 )
            {
                vector Sf( bladePatchNormals[faceI] ); // surface normal
                vector dm( rc - origin ); // moment arm

                vector dFp = rhoRef * Sf * (p.boundaryField()[interfacePatchID][faceI] - pRef) / dr;
                Fp += dFp;
                Mp += dm ^ dFp;

                //vector dFv = mag(Sf) * bladePatchShearStress[faceI] / dr;
                vector dFv = bladePatchShearStress[faceI] / dr;
                Fv += dFv;
                Mv += dm ^ dFv;

                nFacesFound += 1;
            }

        }// end of face loop

        Pstream::gather(Fp, sumOp<vector>()); // These are the DISTRIBUTED loads!
        Pstream::gather(Fv, sumOp<vector>());
        Pstream::gather(Mp, sumOp<vector>());
        Pstream::gather(Mv, sumOp<vector>());

        double Ftot[3], Mtot[3];
        if(Pstream::master())
        {
            for (int i=0; i<3; ++i) 
            {
                Ftot[i] = loadMultiplier*(Fp[i] + Fv[i]);
                Mtot[i] = loadMultiplier*(Mp[i] + Mv[i]);
            }
            //beamDynSetDistributedLoadAtNode(&inode, Ftot, Mtot);
            beamDynSetDistributedLoad(&ig, Ftot, Mtot); // set for each gauss pt
        }

        Pstream::gather(nFacesFound, sumOp<int>());
        Info<< "segment " << ig
            << " with " << nFacesFound << " faces between " << r0 << " " << r1 << ":"
            << " (" << Ftot[0] << " " << Ftot[1] << " " << Ftot[2] << ") " 
            << " (" << Mtot[0] << " " << Mtot[1] << " " << Mtot[2] << ") " 
            << endl;

    } // end loop over beamdyn nodes

