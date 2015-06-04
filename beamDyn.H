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

// We're already using namespace Foam from fvCFD.H
//#ifndef namespaceFoam
//#define namespaceFoam
//    using namespace Foam;
//#endif

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

    scalar  loadMultiplier;

    // global coupling variables

    int nnodes=0;
    int nSurfNodes=0; // from boundary mesh interface patch local points
    
    double currentTime = -1;
    double currentDeltaT = -1;

    vectorList *pos0_ptr,*rot0_ptr; // starting position and orientation of beam nodes
    vectorList *pos_ptr, *rot_ptr;  // current position and orientation of beam nodes
    vectorList *disp_ptr;           // linear displacement at beam nodes
    scalarList *r_ptr;              // spanwise coordinates
    double *h_ptr;                  // shape functions

    // access functions

    vectorList& disp() { return *disp_ptr; };
    vectorList& pos() { return *pos_ptr; };
    vectorList& rot() { return *rot_ptr; };
    scalarList& r() { return *r_ptr; };
    double* h() { return h_ptr; };

    // member functions

    void start( double t0, double dt );
    void stop();
    void update( double t, double dt );

    void updateNodePositions();

    void calculateShapeFunctions( const pointField& pf );

    void updateSectionLoads( const dynamicFvMesh& mesh, 
                             const volScalarField& p, 
                             const incompressible::turbulenceModel& turbulence );

}

#ifdef NoRepository
#   include "beamDyn.C"
#endif

#endif
