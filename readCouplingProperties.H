scalar pi(Foam::constant::mathematical::pi);

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

word interfacePatchName(couplingProperties.lookup("interfacePatch"));

label interfacePatchID =
    mesh.boundaryMesh().findPatchID(interfacePatchName);
Info<< "  Found interfacePatchID for " << interfacePatchName << " : " << interfacePatchID << endl;

if (interfacePatchID < 0)
{
    FatalErrorIn(args.executable())
        << "Problem with finding interface patch"
        << abort(FatalError);
}

// Read reference values

//dimensionedScalar 
//rhoRef( 
//        "rhoRef", dimDensity, 
//        couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0) 
//);
scalar rhoRef( couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0) );
Info<< "  Reference density : " << rhoRef << endl;

label bladeDir( couplingProperties.lookupOrDefault<label>("bladeDir",0) );
Info<< "  Blade axis : " << bladeDir << endl;

//dimensionedVector
//origin( 
//        "origin", dimLength, 
//        couplingProperties.lookupOrDefault<vector>("origin",vector::zero)
//);
vector origin( couplingProperties.lookupOrDefault<vector>("origin",vector::zero) );
Info<< "  Origin (moment reference) : " << origin << endl;

