#include "beamDyn.H"
#include "beamDynInterface.H"

namespace BD
{
    //
    // Routines that directly interface with BeamDyn
    //

    void start( double t0, double dt )
    {
        currentTime = t0;

        if (Pstream::master())
        {
            Info<< "\n================================" << endl;
            Info<< "| Starting BeamDyn" << endl;
            Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" << endl;
            beamDynStart( &t0, &dt );
            //if (canRestart) beamDynRestart( &t0, &dt );
            //else beamDynStart( &t0, &dt );
            beamDynGetNnodes( &nnodes ); // total number of nodes in beam model
            Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

        }

        Pstream::scatter(nnodes);
        r_ptr   = new scalarList(nnodes, 0.0);
        pos0_ptr= new vectorList(nnodes, vector::zero);
        rot0_ptr= new vectorList(nnodes, vector::zero);
        pos_ptr = new vectorList(nnodes, vector::zero);
        rot_ptr = new vectorList(nnodes, vector::zero);
        disp_ptr= new vectorList(nnodes, vector::zero);

        if(Pstream::master())
        {
            if (t0 > 0)
            {
                std::string rstFile("BeamDynState_" + Foam::Time::timeName(t0) + ".dat");
                if (FILE *file = fopen(rstFile.c_str(), "r")) {
                    fclose(file);
                } else {
                    Info<< "Problem opening restart file " << rstFile << endl;
                }   
                   
                // attempt to read data for restart
//                double t, ax, ay, az, bx, by, bz;
//                double ux0, uy0, uz0, tx0, ty0, tz0;
//                double fx0, fy0, fz0, mx0, my0, mz0;
//                int loadSteps=0, dispSteps=0;
//                bool foundLoadStep=false, foundDispStep=false;
//
//                Info<< "Checking for restart data" << endl;
//                while( loadFile >> t >> ax >> ay >> az >> bx >> by >> bz )
//                {
//                    loadSteps += 1;
//                    if (t==t0)
//                    {
//                        fx0 = ax;
//                        fy0 = ay;
//                        fz0 = az;
//                        mx0 = bx;
//                        my0 = by;
//                        mz0 = bz;
//                        foundLoadStep = true;
//                    }
//                }
//                while( dispFile >> t >> ax >> ay >> az >> bx >> by >> bz )
//                {
//                    dispSteps += 1;
//                    if (t==t0)
//                    {
//                        ux0 = ax;
//                        uy0 = ay;
//                        uz0 = az;
//                        tx0 = bx;
//                        ty0 = by;
//                        tz0 = bz;
//                        foundDispStep = true;
//                    }
//                }
//
//                Info<< "Number of load and displacement steps read : "
//                    << loadSteps << " " << dispSteps << endl;
//
//                canRestart = foundLoadStep && foundDispStep && loadSteps==dispSteps;
//                if (canRestart)
//                {
//                    Info<< "Restarting with :" << endl;
//                    Info << "  linear displacement  " << ux0 << " " << uy0 << " " << uz0 << endl;
//                    Info << "  angular displacement " << tx0 << " " << ty0 << " " << tz0 << endl;
//                    Info << "  force distribution  " << fx0 << " " << fy0 << " " << fz0 << endl;
//                    Info << "  moment distribution " << mx0 << " " << my0 << " " << mz0 << endl;
//                }
                beamDynReadState( rstFile.c_str() );

                loadFile.open("load.out", std::ios::in | std::ios::out | std::ios::app);
                dispFile.open("disp.out", std::ios::in | std::ios::out | std::ios::app);
            }
            else
            {
                loadFile.open("load.out", std::ios::out);
                dispFile.open("disp.out", std::ios::out);
            }
            if (!loadFile.is_open()) Info<< "Problem opening load.out???" << endl;
            if (!dispFile.is_open()) Info<< "Problem opening disp.out???" << endl;

            Info<< "Setting precision to " << std::numeric_limits<double>::digits10 << endl;
            loadFile.precision(std::numeric_limits<double>::digits10);
            dispFile.precision(std::numeric_limits<double>::digits10);

        } //if Pstream master

        updateNodePositions();

        for( int inode=0; inode < nnodes; ++inode )
        {
            //for( int i=0; i<3; ++i )
            //{
            //    *pos0_ptr[inode].component(i) = *pos_ptr[inode].component(i);
            //    *rot0_ptr[inode].component(i) = *rot_ptr[inode].component(i);
            //}
            (*pos0_ptr)[inode] = (*pos_ptr)[inode];
            (*rot0_ptr)[inode] = (*rot_ptr)[inode];
        }

        Info<< "BeamDyn initialization complete.\n\n";
    }

    //*********************************************************************************************

    void stop()
    {
        Info<< "================================" << endl;
        Info<< "| Stopping BeamDyn" << endl;
        Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
        if(Pstream::master()) beamDynEnd();
        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << endl;

        delete pos0_ptr;
        delete rot0_ptr;
        delete pos_ptr;
        delete rot_ptr;
        delete disp_ptr;
        delete r_ptr;
        delete [] h_ptr;

        if (Pstream::master)
        {
            loadFile.close();
            dispFile.close();
        }
    }

    //*********************************************************************************************

    void update( double t, double dt )
    {
//        Info<< "================================" << endl;
//        Info<< "| Calling BeamDyn update" << endl;
//        Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
        currentTime = t;
        if(Pstream::master()) 
        {
            beamDynStep( &dt );
        }
//        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << nl << endl;
        updateNodePositions();
    }

    //*********************************************************************************************

    void write( bool writeNow, std::string timeName  )
    {
        if (!writeNow || !Pstream::master()) return;

        std::string fname("BeamDynState_" + timeName + ".dat");
        beamDynWriteState( fname.c_str() );
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Interface calculations
    //

    void updateNodePositions()
    {
        scalarList &r   = *r_ptr;
        vectorList &pos0= *pos0_ptr;
        vectorList &pos = *pos_ptr;
        vectorList &rot0= *rot0_ptr;
        vectorList &rot = *rot_ptr;
        vectorList &disp= *disp_ptr;

//        Info<< "Retrieving node positions for the next iteration" << endl;
        if(Pstream::master())
        {
            dispFile << currentTime;

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

                // get linear/angular displacements
                vector lin_disp( pos[inode] - pos0[inode] );
                vector ang_disp( rot[inode] - rot0[inode] );

                //TODO: handle 3D rotations
                const scalar ang = ang_disp.component(0);   // positive is nose up
                disp[inode].component(0) = 0.0;             // assume no spanwise deformation (in 2D)
                disp[inode].component(2) =                  // chordwise (TE->LE) displacement
                    lin_disp.component(2)*Foam::cos(ang) + lin_disp.component(1)*Foam::sin(ang);
                disp[inode].component(1) =                  // normal displacement
                   -lin_disp.component(2)*Foam::sin(ang) + lin_disp.component(1)*Foam::cos(ang);

                r[inode] = posi[bladeDir];

                dispFile << " " << lin_disp[0] 
                         << " " << lin_disp[1] 
                         << " " << lin_disp[2];
                dispFile << " " << 180/pi*roti[0] 
                         << " " << 180/pi*roti[1] 
                         << " " << 180/pi*roti[2];

            }// loop over beam nodes

            dispFile << std::endl;

        }// if Pstream::master

        Pstream::scatter(r);
        Pstream::scatter(pos); // verified that this works
        Pstream::scatter(rot);
        Pstream::scatter(disp);
    }

    //*********************************************************************************************

    void calculateShapeFunctions( const pointField& pf )
    {
        nSurfNodes = pf.size();
        h_ptr = new double[nSurfNodes*nnodes];
        //double &h = *h_ptr;

        scalarList &r = *r_ptr;
        if( bladeR0 < 0.0 ) bladeR0 = r[0];
        if( bladeR  < 0.0 ) bladeR  = r[nnodes-1];
        //Pout<< "Blade span : " << bladeR0 << " " << bladeR << endl;

        if( nSurfNodes > 0 )
        {
            Pout << "calculating shape functions for " << nSurfNodes << " surface nodes" << endl;

            double s;
            double L_2 = (bladeR-bladeR0)/2.0;
            double hi[nnodes];

            double num, den;
            double GLL[nnodes];
            for( int i=0; i<nnodes; ++i )
            {
                GLL[i] = 2.0*(r[i]-r[0])/(r[nnodes-1]-r[0]) - 1.0;
            }

            forAll( pf, ptI )
            {
                s = ( pf[ptI].component(bladeDir) - bladeR0 ) / L_2 - 1.0;
                //beamDynGetShapeFunctions( &s, hi ); // this only works on the master node...
                //vvvvvvvvvv Code snippet from BeamDyn diffmtc subroutine vvvvvvvvvv
                for( int j=0; j<nnodes; ++j )
                {
                    hi[j] = 0.0;
                    num = 1.0;
                    den = 1.0;
                    if( abs(s-GLL[j]) <= eps )
                    {
                        hi[j] = 1.0;
                    }
                    else
                    {
                        for( int i=0; i<nnodes; ++i )
                        {
                            if( i != j )
                            {
                                den *= (GLL[j] - GLL[i]);
                                num *= (s - GLL[i]);
                            }
                        }
                        hi[j] = num/den;
                    }
                }
                //^^^^^^^^^^^^^^^^^^^ End of code snippet ^^^^^^^^^^^^^^^^^^^^^^^^^^

                for( int inode=0; inode < nnodes; ++inode )
                {
                    h_ptr[ptI*nnodes + inode] = hi[inode];
                }

            }// loop over surface nodes
        }
    }

    //*********************************************************************************************

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
//        Info<< "Calculating surface shear stresses" << endl;
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
        loadFile << currentTime;
//        Info<< "Integrating sectional loads" << endl;
        if(initialInfo) Info<< "Initial info:" << endl;
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

                loadFile << " " << Ftot[0] 
                         << " " << Ftot[1] 
                         << " " << Ftot[2];
                loadFile << " " << Mtot[0] 
                         << " " << Mtot[1] 
                         << " " << Mtot[2];
                
            }

            if (initialInfo)
            {
                Pstream::gather(nFacesFound, sumOp<int>());
                //Info<< "  seg " << ig
                //    << " with " << nFacesFound << " faces btwn " << r0 << " " << r1 << ":"
                //    << " (" << Ftot[0] << " " << Ftot[1] << " " << Ftot[2] << ") " 
                //    << " (" << Mtot[0] << " " << Mtot[1] << " " << Mtot[2] << ") " 
                //    << endl;
                Info<< "segment " << ig
                    << " with " << nFacesFound << " faces"
                    << " between " << r0 << " " << r1 << endl;
            }

        } // end loop over beamdyn nodes

        loadFile << std::endl;

        initialInfo = false;

    }

} // end of BD namespace
