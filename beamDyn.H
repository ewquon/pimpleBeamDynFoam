#ifndef beamDyn_H
#define beamDyn_H

#include "className.H"

// need to include these headers so that compilation of 
// beamDynInterfacePointPatchField doesn't crap out!
#include "fvMesh.H"
#include "fvMatrices.H"
#include "vectorList.H"
#include "dynamicFvMesh.H"
#include "turbulenceModel.H"

#include "constants.H"

#include <fstream>
#include <limits>

namespace BD
{

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

    NamespaceName("BeamDyn");

    const scalar pi(Foam::constant::mathematical::pi);
    const scalar eps(1.0e-8);
    bool first = true;

    std::fstream loadFile;
    std::fstream dispFile;

    // inputs from couplingProperties

    word    interfacePatchName;
    label   interfacePatchID;

    scalar  rhoRef;
    scalar  pRef;
    label   bladeDir;
    scalar  bladeR0;
    scalar  bladeR;
    vector  origin;

    Switch  twoD;

    scalar  loadMultiplier;
 
    vector  prescribed_max_deflection; // only works for y deflections for now
    vector  prescribed_max_rotation;   // only works for x rotations for now, overwrites bending rotation if specified

    // global coupling variables

    int nnodes=0;
    int nSurfNodes=0; // from boundary mesh interface patch local points
    
    double currentTime = -1;
    double currentDeltaT = -1;

    vectorList *pos0_ptr,*rot0_ptr; // starting position and orientation of beam nodes
    vectorList *pos_ptr, *rot_ptr;  // current position and orientation of beam nodes
    vectorList *disp_ptr;           // linear displacement at beam nodes
    vectorList *adisp_ptr;          // angular displacement at beam nodes
    scalarList *r_ptr;              // spanwise coordinates of beam nodes
    double *h_ptr;                  // shape functions [ # CFD nodes * # BD nodes ]

    // access functions

    vectorList& disp() { return *disp_ptr; };
    vectorList& adisp() { return *adisp_ptr; };
    scalarList& r() { return *r_ptr; };
    double* h() { return h_ptr; };

    // member functions

    void start( double t0, double dt );
    void calculateShapeFunctions( const pointField& pf );

    void updateSectionLoads( const dynamicFvMesh& mesh, 
                             const volScalarField& p, 
                             const incompressible::turbulenceModel& turbulence );

    void update( double t, double dt );
    void updatePrescribedDeflection( double t );

    void updateNodePositions();

    void stop();
}

#ifdef NoRepository
#   include "beamDyn.C"
#endif

#endif
