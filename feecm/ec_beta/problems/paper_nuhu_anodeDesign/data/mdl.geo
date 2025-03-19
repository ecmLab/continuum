/// Build a 2D model with micro-electrode

// Parameters
hSE  = 100.0;  // thickness of the SE, in unit um
wSE  = 10.0;  // width of the SE, in unit um
mSE  = 0.5;  // mesh characteristic length for SE out boundaries
mAnd = 0.05;  // mesh characteristic length for micro-probe

//**** Create points **/
// Create points for the boundary
Point(1) = {0, 0, 0, mSE};            Point(2) = {wSE, 0, 0, mSE};
Point(3) = {wSE, hSE, 0, mAnd};       Point(4) = {0, hSE, 0, mAnd}; 

//***** Create lines **/
// Creat lines for the boundary
Line(1)   = {1, 2};     Line(2)   = {2, 3};
Line(3)   = {3, 4};     Line(4)   = {4, 1};

//*** Create loops from lines **/
Curve Loop(5) = {1:4};   // Outer boundary

//*** Create blocks from lines **/
Plane Surface(1) = {5}; // creat block SE

//*** Physical names **/
Physical Curve("SE_ctd")     = {1};
Physical Curve("SE_rgt")     = {2};
Physical Curve("SE_and")     = {3};
Physical Curve("SE_lft")     = {4};
Physical Surface("blockSE")  = {1};
