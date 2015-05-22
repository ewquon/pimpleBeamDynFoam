/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleBeamDynFoam.C

Description
    Transient solver for incompressible, flow of Newtonian fluids
    on a moving mesh using the PIMPLE (merged PISO-SIMPLE) algorithm.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    Interface with NREL FAST/BeamDyn added by Eliot Quon (eliot.quon@nrel.gov).
    Usage:
      - Specify the dynamicMotionSolverFvMesh solver in dynamicMeshDict
      - Apply beamDynInterface[PointPatchVectorField] condition to flexible 
        wall boundaries

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"

// additional includes
#include "beamDyn.H" // BD namespace
#include "scalar.H"
#include "vectorList.H"
#include "pointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
    #include "readTimeControls.H"
    #include "createPcorrTypes.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // additional setup for fsi
    #include "readCouplingProperties.H"
    double t0 = runTime.startTime().value();
    double dt = runTime.deltaT().value();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    BD::start( t0, dt );

    #include "beamDynVars.H" // TODO: remove this, use BD namespace instead!

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "CourantNo.H"

        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Displacements are updated through the beamDynInterface boundary condition
        // Note: there should not be any displacement for the first time step
        mesh.update();

        // Calculate absolute flux from the mapped surface velocity
        phi = mesh.Sf() & Uf;

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        // additional fsi steps

        Info<< "\nCalculating sectional loads for BeamDyn" << endl;
        BD::updateSectionLoads( mesh, p, turbulence );

        BD::update( runTime.deltaT().value() );

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    BD::stop();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
