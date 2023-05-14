/// Build a 2D model with micro-electrode

SetFactory("OpenCASCADE");

// Parameters
SEx  = 5000;  // length of the SE, in unit um
SEy  = 500;   // height of the SE, in unit um
andX = 10;   // height of the SE, in unit um
eldH = 50;  // height of the electrodes, fixed to 5um
mSE  = 10;  // mesh characteristic length for SE out boundaries
mCtd = 5;  // mesh characteristic length for cathode
mAnd = 2;  // mesh characteristic length for cathode

//**** Create the Electrodes **/
Rectangle(1) = {-andX/2,SEy,0,andX,eldH};
Rectangle(2) = {-SEx/2, -eldH,0,SEx,eldH};
//**** Create the SE **/
Rectangle(3) = {-SEx/2,0,0,SEx,SEy};

//*** Coherence interface **/
v() = BooleanFragments { Surface{1}; Delete; }{ Surface{3}; Delete; };
v() = BooleanFragments { Surface{2}; Delete; }{ Surface{3}; Delete; };

//*** Mesh control **/
Characteristic Length{ PointsOf{ Curve{1,3}; } } = mAnd;
Characteristic Length{ PointsOf{ Curve{5}; } } = mCtd;
Characteristic Length{ PointsOf{ Curve{7,9,12}; } } = mSE;

//*** Physical names **/
Physical Curve("And_bottom") = {1};
Physical Curve("And_right")  = {2};
Physical Curve("And_top")    = {3};
Physical Curve("And_left")   = {4};
Physical Curve("Ctd_bottom") = {5};
Physical Curve("Ctd_right")  = {6};
Physical Curve("Ctd_top")    = {7};
Physical Curve("Ctd_left")   = {8};
Physical Curve("SE_right")   = {9};
Physical Curve("SE_top")     = {10,11};
Physical Curve("SE_left")    = {12};
Physical Surface("blockAnd") = {1};
Physical Surface("blockCtd") = {2};
Physical Surface("blockSE") = {3};
