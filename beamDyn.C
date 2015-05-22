#include "beamDyn.H"
#include "beamDynInterface.H"

namespace BD
{

    //
    // Routines that directly interface with BeamDyn
    //

    void start( double t0, double dt )
    {
        Info<< "\n================================" << endl;
        Info<< "| Starting BeamDyn" << endl;
        Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" << endl;
        if(Pstream::master())
        {
            beamDynStart( &t0, &dt );
            beamDynGetNnodes( &nnodes ); // total number of nodes in beam model
        }
        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

        Pstream::scatter(nnodes);
        pos_ptr = new vectorList(nnodes, vector::zero);
        rot_ptr = new vectorList(nnodes, vector::zero);
        r_ptr   = new scalarList(nnodes, 0.0);

        updateNodePositions();
    }

    void stop()
    {
        Info<< "================================" << endl;
        Info<< "| Stopping BeamDyn" << endl;
        Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
        if(Pstream::master()) beamDynEnd();
        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << endl;
        delete pos_ptr;
        delete rot_ptr;
        delete r_ptr;
    }

    void update( double dt )
    {
        Info<< "\n================================" << endl;
        Info<< "| Calling BeamDyn update" << endl;
        Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
        if(Pstream::master()) 
        {
            beamDynStep( &dt );
        }
        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << endl;
        updateNodePositions();
    }

    //
    // Interface calculations
    //

    void updateNodePositions()
    {
//        Foam::vectorList pos( nnodes, Foam::vector::zero );    // current position (for calling beamDynGetNodePosition)
//        Foam::vectorList rot( nnodes, Foam::vector::zero );    // current orientation (for calling beamDynGetNodePosition)
//        Foam::scalarList r( nnodes, 0.0 );                     // spanwise coordinate, used in updateSectionLoads.H
        vectorList &pos = *pos_ptr;
        vectorList &rot = *rot_ptr;
        scalarList &r   = *r_ptr;

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
    }

    void updateSectionLoads( const dynamicFvMesh& mesh, 
                             const volScalarField& p, 
                             const incompressible::turbulenceModel& turbulence )
    {
        scalarList &r = *r_ptr;
        scalar p0( pRef / rhoRef );

        // setup arrays, pointers
        double r0, r1;
        const polyPatch& bladePatch = mesh.boundaryMesh()[interfacePatchID];
        const vectorField& bladePatchNormals = mesh.Sf().boundaryField()[interfacePatchID];

        // calculate shear stress
        Info<< "Calculating surface shear stresses" << endl;
        //const volSymmTensorField Reff(turbulence->devReff());
        const volSymmTensorField Reff(turbulence.devReff());
        vectorField bladePatchShearStress = 
            (
                -mesh.Sf().boundaryField()[interfacePatchID]
                /mesh.magSf().boundaryField()[interfacePatchID]
            ) & Reff.boundaryField()[interfacePatchID];

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
                    vector dm( rc - origin );

                    vector dFp = rhoRef * Sf * (p.boundaryField()[interfacePatchID][faceI] - p0) / dr;
                    Fp += dFp;
                    Mp += dm ^ dFp;

                    vector dFv = mag(Sf) * bladePatchShearStress[faceI] / dr;
                    Fv += dFv;
                    Mv += dm ^ dFv;

                    nFacesFound += 1;
                }

            }// end of face loop
            //Pout << "nFacesFound= " << nFacesFound << " between " << r0 << " " << r1 << endl;

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
                beamDynSetDistributedLoad(&ig, Ftot, Mtot);
            }

            Pstream::gather(nFacesFound, sumOp<int>());
            Info<< "segment " << ig
                << " with " << nFacesFound << " faces between " << r0 << " " << r1 << ":"
                << " (" << Ftot[0] << " " << Ftot[1] << " " << Ftot[2] << ") " 
                << " (" << Mtot[0] << " " << Mtot[1] << " " << Mtot[2] << ") " 
                << endl;

        } // end loop over beamdyn nodes
    }

}
