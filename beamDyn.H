#ifndef beamDyn_H
#define beamDyn_H

#include "className.H"
#include "mathematicalConstants.H"
#include "vectorList.H"

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
    
    vectorList *pos_ptr, *rot_ptr;  // current position and orientation of beam nodes
    scalarList *r_ptr;              // spanwise coordinates

    // access functions

    vectorList& pos() { return *pos_ptr; };
    vectorList& rot() { return *rot_ptr; };
    scalarList& r() { return *r_ptr; };

    // member functions

    void start( double t0, double dt );
    void stop();
    void update( double dt );

    void updateNodePositions();

    void updateSectionLoads();

}

#ifdef NoRepository
#   include "beamDyn.C"
#endif

#endif
