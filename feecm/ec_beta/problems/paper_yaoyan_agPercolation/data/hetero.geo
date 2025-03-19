lx   = 10;       // width of the BL, in unit um
ly   = 20;       // thickness of the BL, in unit um
dia  = 0.5;       // diameter of the pore located at the middle, in unit um

m0BL  = 0.1;     // mesh characteristic length for the bl miec
m1Pr  = 0.01;    // mesh characteristic length for Ag particles or pores

//**** Create BL geometry**/
Point(1)    = {-lx/2, -ly/2, 0, m0BL};
Point(2)    = { lx/2, -ly/2, 0, m0BL};
Point(3)    = { lx/2,  ly/2, 0, m0BL};
Point(4)    = {-lx/2,  ly/2, 0, m0BL};
Line(1)     = {1,2};
Line(2)     = {2,3};
Line(3)     = {3,4};
Line(4)     = {4,1};
Curve Loop(51) = {1,2,3,4};

//*** Create pore **/
Point(5)    = {-dia/2-lx/4, 0,       0, m1Pr};
Point(6)    = { -lx/4,      -dia/2,  0, m1Pr};
Point(7)    = { dia/2-lx/4, 0,       0, m1Pr};
Point(8)    = { -lx/4,      dia/2,   0, m1Pr};
Point(9)    = { -lx/4,      0,       0, m1Pr};
Circle(5)   = {5,9,6};
Circle(6)   = {6,9,7};
Circle(7)   = {7,9,8};
Circle(8)   = {8,9,5};
Curve Loop(52) = {5,6,7,8};

//*** Create Ag particle **/
Point(10)   = {-dia/2+lx/4, 0,       0, m1Pr};
Point(11)   = { lx/4,       -dia/2,  0, m1Pr};
Point(12)   = { dia/2+lx/4, 0,       0, m1Pr};
Point(13)   = { lx/4,       dia/2,   0, m1Pr};
Point(14)   = { lx/4,       0,       0, m1Pr};
Circle(9)   = {10,14,11};
Circle(10)  = {11,14,12};
Circle(11)  = {12,14,13};
Circle(12)  = {13,14,10};
Curve Loop(53) = {9,10,11,12};

//*** Create blocks for the BL and pore **/
Plane Surface(1) = {51,52,53};

//*** Mesh control
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockBL_btm") = {1};
//+
Physical Curve("blockBL_rgt") = {2};
//+
Physical Curve("blockBL_top") = {3};
//+
Physical Curve("blockBL_lft") = {4};
//+
Physical Curve("pore")     = {5,6,7,8};
//+
Physical Curve("Ag")       = {9,10,11,12};
//+
Physical Surface("blockBL") = {1};
//+
//Physical Surface("blockAg") = {2};
