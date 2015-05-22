//scalar pi(Foam::constant::mathematical::pi);

Info << "\nReading coupling properties" << endl;
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Read interface patch info
// TODO: only one interface patch for now

//word interfacePatchName( couplingProperties.lookup("interfacePatch") );
//BD::interfacePatchName = couplingProperties.lookup("interfacePatch");
couplingProperties.lookup("interfacePatch") >> BD::interfacePatchName;
BD::interfacePatchID = mesh.boundaryMesh().findPatchID(BD::interfacePatchName);
Info<< "  Found interfacePatchID for " << BD::interfacePatchName 
    << " : " << BD::interfacePatchID << endl;

if (BD::interfacePatchID < 0)
{
    FatalErrorIn(args.executable())
        << "Problem with finding interface patch"
        << abort(FatalError);
}

const polyPatch& interfacePatch = mesh.boundaryMesh()[BD::interfacePatchID];

// Read reference values

//dimensionedScalar 
//rhoRef( 
//        "rhoRef", dimDensity, 
//        couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0) 
//);
//scalar rhoRef( couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0) );
BD::rhoRef = couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0);
BD::pRef = pRefValue;
Info<< "  Reference density  : " << BD::rhoRef << endl;
Info<< "  Reference pressure : " << BD::pRef << endl;

//label bladeDir( couplingProperties.lookupOrDefault<label>("bladeDir",0) );
BD::bladeDir = couplingProperties.lookupOrDefault<label>("bladeDir",0);
Info<< "  Blade axis : " << BD::bladeDir << endl;

//scalar bladeR( couplingProperties.lookupOrDefault("R",-1) );
//scalar bladeR0( couplingProperties.lookupOrDefault("R0",-1) );
BD::bladeR0 = couplingProperties.lookupOrDefault("R0",-1);
BD::bladeR  = couplingProperties.lookupOrDefault("R",-1);

//dimensionedVector
//origin( 
//        "origin", dimLength, 
//        couplingProperties.lookupOrDefault<vector>("origin",vector::zero)
//);
//vector origin( couplingProperties.lookupOrDefault<vector>("origin",vector::zero) );
BD::origin = couplingProperties.lookupOrDefault<vector>("origin",vector::zero);
Info<< "  Origin (moment reference) : " << BD::origin << endl;

//scalar loadMultiplier( couplingProperties.lookupOrDefault<scalar>("loadMultiplier",1.0) );
BD::loadMultiplier = couplingProperties.lookupOrDefault<scalar>("loadMultiplier",1.0);
Info<< "  Load multiplier (FOR DEVELOPMENT) : " << BD::loadMultiplier << endl;

