/// Build a 2D model with micro-electrode

//SetFactory("OpenCASCADE");

// Parameters
hSE  = 20.0;  // thickness of the SE, in unit um
wSE  = 8.0;  // width of the SE, in unit um
hLi  = 4.4006;    // initial height of the Li in the pore, in unit um
rp   = 0.05;  // the round radius at corner point, in unit um
mSE  = 0.5;  // mesh characteristic length for SE out boundaries
mAnd = 0.05;  // mesh characteristic length for micro-probe

//**** Create points **/
// Create points for the boundary
Point(1) = {0, 0, 0, mSE};            Point(2) = {wSE, 0, 0, mSE};
Point(3) = {wSE, hSE-hLi, 0, mAnd};   Point(4) = {wSE, hSE-rp, 0, mAnd}; 
Point(5) = {wSE-rp, hSE-rp, 0, mAnd}; Point(6) = {wSE-rp, hSE, 0, mAnd}; 
Point(7) = {0, hSE, 0, mAnd};

//***** Create lines **/
// Creat lines for the boundary
Line(1)   = {1, 2};     Line(2)   = {2, 3};
Line(3)   = {3, 4};     Circle(4) = {4, 5, 6};  
Line(5)   = {6, 7};     Line(6)   = {7, 1};

//*** Create loops from lines **/
Curve Loop(6) = {1, 2, 3, 4, 5, 6};   // Outer boundary

//*** Create blocks from lines **/
Plane Surface(1) = {6}; // creat block SE

//*** Physical names **/
Physical Curve("SE_btm")     = {1};
Physical Curve("SE_right")   = {2};
//Physical Curve("SE_top")     = {3};
Physical Curve("SE_left")    = {6};
Physical Curve("And_btm")    = {3, 4, 5};
Physical Surface("blockSE")  = {1};
