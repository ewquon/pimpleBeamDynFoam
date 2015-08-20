#include "beamDyn.H"
#include "beamDynInterface.H"

namespace BD
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Routines that directly interface with BeamDyn
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    void start( double t0, double dt )
    {
        currentTime = t0;
        currentDeltaT = dt;

        if (Pstream::master())
        {
            Info<< "\n================================" << endl;
            Info<< "| Starting BeamDyn" << endl;
            Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" << endl;
            beamDynStart( &t0, &dt );
            beamDynGetNnodes( &nnodes ); // total number of nodes in beam model
            Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
        }

        // initialize arrays for storing configuration
        Pstream::scatter(nnodes);
        r_ptr    = new scalarList(nnodes, 0.0);
        pos0_ptr = new vectorList(nnodes, vector::zero);
        rot0_ptr = new vectorList(nnodes, vector::zero);
        pos_ptr  = new vectorList(nnodes, vector::zero);
        rot_ptr  = new vectorList(nnodes, vector::zero);
        disp_ptr = new vectorList(nnodes, vector::zero);
        adisp_ptr= new vectorList(nnodes, vector::zero);

        // perform restart read of saved state data if necessary
        // open disp.out and load.out for writing later
        if(Pstream::master())
        {
            if (t0 > 0)
            {
                Info<< "Time = " << t0 << ", attempting restart" << endl;
                std::string rstFile("BeamDynState_" + Foam::Time::timeName(t0) + ".dat");
                if (FILE *file = fopen(rstFile.c_str(), "r")) {
                    fclose(file);
                } else {
                    Info<< "Problem opening restart file " << rstFile << endl;
                }   
                   
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

            //Info<< "Setting precision to " << std::numeric_limits<double>::digits10 << endl;
            //loadFile.precision(std::numeric_limits<double>::digits10);
            //dispFile.precision(std::numeric_limits<double>::digits10);
            loadFile.precision(8);
            dispFile.precision(8);

            // get initial configuration
            double posi[3], roti[3];
            Info<< "Initial linear/angular position:" << endl;
            for( int inode=0; inode < nnodes; ++inode )
            {
                beamDynGetInitNode0Position( &inode, posi, roti );
                for(int i=0; i < 3; ++i) {
                    (*pos0_ptr)[inode][i] = posi[i];
                    (*rot0_ptr)[inode][i] = roti[i];
                }
                Info<< " " << posi[0] 
                    << " " << posi[1] 
                    << " " << posi[2]
                    << " " << 180/pi*roti[0] 
                    << " " << 180/pi*roti[1] 
                    << " " << 180/pi*roti[2]
                    << endl;
            }

        } //if Pstream master

        updateNodePositions(); // this should write out either the initial configuration (0's)
                               // or the restart configuration

//        // save initial configuration, used by updateNodePositions() in subsequent iterations
//        for( int inode=0; inode < nnodes; ++inode )
//        {
//            (*pos0_ptr)[inode] = (*pos_ptr)[inode];
//            (*rot0_ptr)[inode] = (*rot_ptr)[inode];
//        }

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
        delete adisp_ptr;
        delete r_ptr;
        delete [] h_ptr;

        if (Pstream::master())
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
        currentDeltaT = dt;
        if(Pstream::master()) 
        {
            beamDynStep( &dt );
        }
//        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << nl << endl;
        updateNodePositions();
    }
    //*********************************************************************************************

    void updatePrescribedDeflection( double t )
    {
        if(Pstream::master()) 
        {
            // prescribed deflection for testing
            // y = ymax * (x/L)^2
            //
            // - at t=1.0, deflection is 100%
            // - calculate shortening of the blade, normalized by blade length:
            // b = (3*a**4 + sqrt(9*a**8+2*a**6))**(1./3.)
            // xmax = b/2**(2./3.)/a**2 - 1/2**(1./3.)/b

            double L = bladeR - bladeR0;
            double a = 2*prescribed_max_deflection[1]*t; // prescribed deflection only in y-dir for now
            double ang = prescribed_max_rotation[0]*currentDeltaT; // rotation increment
            Info<< "current delta t : " << currentDeltaT << endl;
            if (a != 0)
            {
                double b = Foam::pow( 3*Foam::pow(a,4) 
                                    + Foam::sqrt( 9*Foam::pow(a,8) + 2*Foam::pow(a,6) )
                                 ,1.0/3.0);
                double xmax = L * ( b/Foam::pow(2,2.0/3.0)/pow(a,2) - 1.0/(b*pow(2,1.0/3.0)) );
                double ymax = L * prescribed_max_deflection[1]*t;

                double x[3], tmp[3], u[3];//, p[3];
                for (int i=0; i<3; ++i)
                {
                    u[i] = 0;
                }

                Info<< "Prescribed tip position at time " << t 
                    << " : x/L,y/L = " << xmax/L << " " << ymax/L 
                    << " with twist increment = " << ang
                    << endl;

                for (int inode=0; inode<nnodes; ++inode)
                {
                    //beamDynGetInitNode0Position( &inode, x0, tmp );
                    beamDynGetNode0Position( &inode, x, tmp );
                    //double xi = x[bladeDir]/L;
                    double s = (*pos0_ptr)[inode][bladeDir] / L; //parametric coordinate from 0 to 1

                    // NOTE: prescribed deflection is only in y-dir for now
                    //       bladeDir is implied to be x-dir
                    u[bladeDir] = (xmax-L)*s; //xi * xmax - (*pos0_ptr)[inode][bladeDir];
                    u[1] = ymax*s*s - (*pos0_ptr)[inode][1];
                    u[2] = 0.0;

    // DEBUG
    //                Info<< "  initial position is " << (*pos0_ptr)[inode] 
    //                    << "  (s=" << s << ")"
    //                    << endl;
    //
    //                Info<< "  setting disp for node " << inode+1
    //                    << " at " << x[0] << " " << x[1] << " " << x[2]
    //                    << " to " << u[0] << " " << u[1] << " " << u[2]
    //                    << endl;

                    beamDynSetDisplacement( &inode, u );

                    // Since we're prescribing motion and not actually running the BeamDyn solver,
                    // the rotation matrix is never actually updated. We manually set it here.
                    //double dydx = a*xi;
                    //double ang = Foam::atan(a*xi);
                    ang = Foam::atan(a*s);
                    beamDynSetZRotationMatrix( &inode, &ang );

                }
            }

            // For now, just overwrite the rotation matrix if we're doing a test with x-axis rotation
            if (ang != 0)
            {
                Info<< "Specifying twist increment = " << ang << endl;
                for (int inode=0; inode<nnodes; ++inode)
                {
//                    ang = prescribed_max_rotation[0]*t + (*rot0_ptr)[inode][0]; // this is a rotation angle, not the absolute orientation!
//                    ang = prescribed_max_rotation[0]*currentDeltaT + (*rot0_ptr)[inode][0];
//                    Info<< "Prescribing twist angle " << ang*180/pi << endl;
                    beamDynSetXRotationMatrix( &inode, &ang );
                }
            }
        }

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
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Retrieve disp array from the BeamDyn library
    // TODO: clean this up?
    // - updates pos/rot, used to calculate disp
    // - updates disp, accessed through BD::disp() in beamDynInterfacePointPatch::updateCoeffs()
    // - updates r, used by updateSectionLoads()
    // - writes displacement at the starting time step to disp.out
    void updateNodePositions()
    {
        Info<< "Updating node positions from BeamDyn" << endl;

        scalarList &r    = *r_ptr;
        vectorList &pos0 = *pos0_ptr;
        //vectorList &rot0 = *rot0_ptr; // not used
        vectorList &pos  = *pos_ptr;
        vectorList &rot  = *rot_ptr;
        vectorList &disp = *disp_ptr;
        vectorList &adisp= *adisp_ptr;

//        Info<< "Retrieving node positions for the next iteration" << endl;
        if(Pstream::master())
        {
            if (first) Info<< "Initial displacements from BeamDyn lib [m,deg]: " << endl;
            else dispFile << currentTime;

            // --loop over nodes in the BeamDyn blade model (assumed single element)
            //   TODO: handle multiple elements
//            double posi[3], roti[3];
            double lin_disp[3], ang_disp[3];
            //double ang;
//            double R[9]; //rotation matrix from BeamDyn
            for( int inode=0; inode<nnodes; ++inode ) 
            {

                // get node linear/angular displacement [m, rad]
                //   returns rotation_angle - initial_rotation_angle
                //   initial_rotation_angle should be stored in (*rot0_ptr) during start()
                //   => this returns the wrong angular displacement at the zeroth step
                //      but is this ever used???
                beamDynGetNode0Displacement( &inode, lin_disp, ang_disp );
// DEBUG
//                Info<< "node " << inode << " : lin_disp " 
//                    << "(" << lin_disp[0] << " " << lin_disp[1] << " " << lin_disp[2] << ")" 
//                    << endl;

// 2D operation, e.g. wingMotion case
//                ang = ang_disp[0];   // positive is nose up
//                disp[inode].component(0) = 0.0;                                                      // assume no spanwise deformation (in 2D)
//                disp[inode].component(2) =  lin_disp[2]*Foam::cos(ang) + lin_disp[1]*Foam::sin(ang); // chordwise (TE->LE) displacement
//                disp[inode].component(1) = -lin_disp[2]*Foam::sin(ang) + lin_disp[1]*Foam::cos(ang); // normal displacement
//                r[inode] = pos0[inode][bladeDir] + disp[inode].component(bladeDir);

// ROTATION SHOULD BE APPLIED FOR THE POINTS ON THE INTERFACE PATCH!!!
//                beamDynGetNode0RotationMatrix( &inode, R );
// DEBUG
// these are identity before the first iteration and approximately the identity matrix afterwards
//                Info << "DEBUG R matrix " << inode << " : \n";
//                Info << " R=[" << R[0] << " " << R[1] << " " << R[2] << "; ...\n";
//                Info << "    " << R[3] << " " << R[4] << " " << R[5] << "; ...\n";
//                Info << "    " << R[6] << " " << R[7] << " " << R[8] << "]\n";
//
//                disp[inode] = vector::zero;
//                for( int j=0; j<3; ++j ) {
//                    for( int i=0; i<3; ++i )
//                    {
//                        disp[inode].component(j) += R[3*j+i] * lin_disp[i];
//                    }
//
//                    //disp[inode].component(j) = lin_disp[j];
//                }
//                if(twoD) disp[inode].component(bladeDir) = 0.0;
                for( int i=0; i<3; ++i )
                {
                    disp[inode].component(i) = lin_disp[i];
                    adisp[inode].component(i)= ang_disp[i];
                }

                // used for sectional loads calculation
                r[inode] = pos0[inode][bladeDir] + disp[inode].component(bladeDir);

// DEBUG
//                Info<< "DEBUG node " << inode << " (lin_disp) :" 
//                    << " " << lin_disp[0]
//                    << " " << lin_disp[1]
//                    << " " << lin_disp[2]
//                    << endl;
//                Info<< "DEBUG node " << inode << " (OLD) :" 
//                    << " " << 0.0
//                    << " " << -lin_disp[2]*Foam::sin(ang_disp[0]) + lin_disp[1]*Foam::cos(ang_disp[0]) // normal displacement
//                    << " " <<  lin_disp[2]*Foam::cos(ang_disp[0]) + lin_disp[1]*Foam::sin(ang_disp[0]) // chordwise (TE->LE) displacement
//                    << endl;
//                Info<< "DEBUG node " << inode << " (NEW) :" 
//                    << " " << disp[inode].component(0)
//                    << " " << disp[inode].component(1)
//                    << " " << disp[inode].component(2)
//                    << endl;


                if (first) // print out initial displaced config, either 0's or (hopefully) repeated on restart
                {
                    //NOTE: for prescribed motion, since the beamdyn solve routine isn't called, the rotation parameters
                    //  are not properly updated so ang_disp will be inaccurate; just use the rotation matrix R instead.
                    Info<< " " << lin_disp[0] 
                        << " " << lin_disp[1] 
                        << " " << lin_disp[2]
                        << " " << 180/pi*ang_disp[0] 
                        << " " << 180/pi*ang_disp[1] 
                        << " " << 180/pi*ang_disp[2] << endl;
                }
                else // write subsequent displacements to file
                {
                    dispFile << " " << lin_disp[0] 
                             << " " << lin_disp[1] 
                             << " " << lin_disp[2];
                    dispFile << " " << 180/pi*ang_disp[0] 
                             << " " << 180/pi*ang_disp[1] 
                             << " " << 180/pi*ang_disp[2];
                }

            }// loop over beam nodes

            if (first) Info<< endl;
            else dispFile << std::endl;

        }// if Pstream::master

        Pstream::scatter(r);
        Pstream::scatter(pos); // verified that this works
        Pstream::scatter(rot);
        Pstream::scatter(disp);
        Pstream::scatter(adisp);
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
        Info<< "r: " << r << endl;
        //Pout<< "Blade span : " << bladeR0 << " " << bladeR << endl;

        if( nSurfNodes > 0 )
        {
            Pout << "calculating shape functions for " << nSurfNodes << " surface nodes" << endl;

            // DEBUG
            scalar hmin( 9e9);
            scalar hmax(-9e9);
            scalar smin( 9e9);
            scalar smax(-9e9);

            double s;
            double L_2 = (bladeR-bladeR0)/2.0;
            double hi[nnodes];

            double num, den;
            double GLL[nnodes];
            Info<< "GLL pts : "; // DEBUG
            for( int i=0; i<nnodes; ++i )
            {
                GLL[i] = 2.0*(r[i]-r[0])/(r[nnodes-1]-r[0]) - 1.0;
                Info<< " " << GLL[i];
            }
            Info<< endl;

            forAll( pf, ptI )
            {
                s = ( pf[ptI].component(bladeDir) - bladeR0 ) / L_2 - 1.0;
                smin = min(smin,s);
                smax = max(smax,s);
                //beamDynGetShapeFunctions( &s, hi ); // this only works on the master node...
                //vvvvvvvvvv Code snippet from BeamDyn diffmtc subroutine vvvvvvvvvv
                for( int j=0; j<nnodes; ++j )
                {
                    hi[j] = 0.0;
                    num = 1.0;
                    den = 1.0;
                    if( fabs(s-GLL[j]) <= eps )
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

                scalar hsum(0);
                for( int inode=0; inode < nnodes; ++inode )
                {
                    h_ptr[ptI*nnodes + inode] = hi[inode];
                    hsum += hi[inode];
                }
                //Info<< hsum << endl; between 3 and 5?!?
                hmin = min(hmin,hsum);
                hmax = max(hmax,hsum);

            }// loop over surface nodes
//            Pout<< "-- sum(h) : [ " << hmin << ", " << hmax << " ]" 
//                << " for s : [ " << smin << ", " << smax << " ]"
//                << endl;
        }
    }

    //*********************************************************************************************

    void updateSectionLoads( const dynamicFvMesh& mesh, 
                             const volScalarField& p, 
                             const incompressible::turbulenceModel& turbulence )
    {
        Info<< "Calculating section loads for BeamDyn" << endl;

        scalarList &r = *r_ptr;
        scalar p0( pRef / rhoRef );

        // setup arrays, pointers
        double r0, r1;
        const polyPatch& bladePatch = mesh.boundaryMesh()[interfacePatchID];
        const vectorField& bladePatchNormals = mesh.Sf().boundaryField()[interfacePatchID];

        // calculate shear stress
        //   note: devReff returns the effective stress tensor including the laminar stress
        //   note: face normals point _outside_ the computational domain

//        Info<< "Calculating surface shear stresses" << endl;
//        const volSymmTensorField Reff(turbulence.devReff());
//        vectorField bladePatchShearStress = 
//            (
//                -mesh.Sf().boundaryField()[interfacePatchID]
//                /mesh.magSf().boundaryField()[interfacePatchID]
//            ) & Reff.boundaryField()[interfacePatchID];

        // Face normals point into solid surface, i.e., outward from fluid volume, 
        // i.e. the direction the fluid is pushing on the wall.
        // This matches the 'forces' function object implementation 
        // in Foam::forces::calcForcesMoment() at
        //   ~/OpenFOAM/OpenFOAM-2.3.1/src/postProcessing/functionObjects/forces/forces/forces.C
        // - also, no need to normalize by magSf since we multiply by mag(Sf) later
        const volSymmTensorField Reff(turbulence.devReff());
        //THIS DOESN'T COMPILE WITH CLANG:
        //vectorField bladePatchShearStress = 
        //    mesh.Sf().boundaryField()[interfacePatchID]
        //    & Reff.boundaryField()[interfacePatchID];
        vectorField bladePatchShearStress( bladePatchNormals & Reff.boundaryField()[interfacePatchID] );

        //
        // --loop over nodes in the BeamDyn blade model, assumed single element
        //   i.e., nnodes = nodes_elem = order_elem+1 = ngp+1
        //
        //loadFile << currentTime; // at this point, still equal to t at beginning of time step
        loadFile << currentTime + currentDeltaT;
//        Info<< "Integrating sectional loads" << endl;
        if(first) Info<< "Initial info:" << endl;
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

                // Update F_foam and M_foam in BeamDyn library
                //beamDynSetDistributedLoadAtNode(&inode, Ftot, Mtot);
                beamDynSetDistributedLoad(&ig, Ftot, Mtot);

                loadFile << " " << Ftot[0] 
                         << " " << Ftot[1] 
                         << " " << Ftot[2];
                loadFile << " " << Mtot[0] 
                         << " " << Mtot[1] 
                         << " " << Mtot[2];
            }

            if (first)
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

        first = false;

    }

} // end of BD namespace
