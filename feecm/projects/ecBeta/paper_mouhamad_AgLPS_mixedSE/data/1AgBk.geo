lx   = 16;       // width of the SE, in unit um
ly   = 20;       // thickness of the SE, in unit um
dh   = 10;       // length of the pre-existing crack, in unit um
dw   = 2;        // width of the pre-existing crack located at the bottom middle, in unit um
hLi  = 5;        // length of the current Li metal, in unit um
hAg  = hLi/2;    // location of the Ag, in unit um
dAg  = 0.1;      // size of the Ag, in unit um

m0SE  = 0.2;    // mesh characteristic length for corner points of SE
m1SE  = 0.002;   // mesh characteristic length for defect points of SE

//**** Create SE geometry**/
Point(1)    = {dw/2, 0, 0, m1SE};
Point(2)    = {lx/2, 0,  0, m0SE};
Point(3)    = {lx/2, ly, 0, m0SE};
Point(4)    = {0,    ly, 0, m0SE};
Point(5)    = {0,    dh, 0, m1SE};
Point(6)    = {0,    dh-dw/2, 0, m1SE};
Point(7)    = {dw/2, dh-dw/2, 0, m1SE};
Point(8)    = {dw/2, hLi, 0, m1SE};
Point(9)    = {dw/2, hAg+dAg/2, 0, m1SE};
Point(10)   = {dw/2, hAg-dAg/2, 0, m1SE};
//Point(11)   = {dw/2, hAg, 0, m1SE};
Line(1)     = {1,2};
Line(2)     = {2,3};
Line(3)     = {3,4};
Line(4)     = {4,5};
Circle(5)   = {5,6,7};
Line(6)     = {7,8};
Line(7)     = {8,9};
//Circle(8)   = {10,11,9};
Line(8)     = {9,10};
Line(9)     = {10,1};

//*** Create loops for SE **/
//Curve Loop(10) = {1:7,-8,9};
Curve Loop(10) = {1:9};

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {10};

//*** Mesh control
//Mesh.RecombineAll = 1;
//Mesh.Algorithm = 6; // Frontal-Delaunay
//Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
//Physical Curve("blockSE_btm")   = {1};
//+
Physical Curve("blockSE_rgt") = {2};
//+
Physical Curve("blockSE_top")   = {3};
//+
Physical Curve("blockSE_lft")  = {4};
//+
Physical Curve("liMetal")   = {1,7,9};
//+
Physical Curve("ag")   = {8};
//+
Physical Surface("blockSE") = {1};

