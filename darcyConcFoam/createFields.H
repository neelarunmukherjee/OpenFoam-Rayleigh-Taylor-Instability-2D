	Info<< "Reading field p\n" << endl;

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
		dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
    );

	surfaceScalarField phi ("phi", linearInterpolate(U) & mesh.Sf());

    Info<< "Reading field c\n" << endl;

    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity Dc\n" << endl;

    dimensionedScalar Dc
    (
        transportProperties.lookup("Dc")
    );
	
	Info<< "Reading diffusivity c0\n" << endl;
	
	dimensionedScalar c0
    (
        transportProperties.lookup("c0")
    );
	
	Info<< "Reading diffusivity del_rho\n" << endl;
	
		dimensionedScalar del_rho
    (
        transportProperties.lookup("del_rho")
    );
 
	Info<< "Reading diffusivity alpha\n" << endl;
	
	dimensionedScalar a
    (
     transportProperties.lookup("a")
    );
	
	Info<< "Reading diffusivity S\n" << endl;
	
	/*dimensionedScalar S
    (
     transportProperties.lookup("S")
    );*/
    
	Info<< "Reading diffusivity a small number\n" << endl;
	
    dimensionedScalar delS
    (
     transportProperties.lookup("delS")
    );
	
	Info<< "Reading diffusivity k\n" << endl;
	
    dimensionedScalar k
    (
        transportProperties.lookup("k")
    );

	Info<< "Reading diffusivity mu\n" << endl;

    dimensionedScalar mu
    (
        transportProperties.lookup("mu")
    );


