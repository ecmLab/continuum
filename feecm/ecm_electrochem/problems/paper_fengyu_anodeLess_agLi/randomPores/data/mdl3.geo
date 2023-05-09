/// Model 1: Build a model with three void vertically distributed
///  With the first hole filled in with a thin layer of Li metal

// Parameters 
lx   = 1;  // length of the model, in unit um
ly   = 1;  //height of the model, in unit um
cntR = 0.05; // Radius of the hole
dHole = ly/4; // Distance of each hole
cntx = lx/2; // hole always stay at the middle
cnty1 = ly - dHole; // y coordinate of the first hole center
cnty2 = ly - 2*dHole; // y coordinate of the second hole center
cnty3 = ly - 3*dHole; // y coordinate of the second hole center

m1    = 0.01; // mesh characteristic length for boundary points
m2    = 0.001; // mesh characteristic length for interface points
m3    = 0.0001; // mesh characteristic length for inner points
//thk   = 0.001; // Thickness of the deposited Li metal

//**** Create points **/
// Create points for the boundary
Point(1) = {0, 0, 0, m1};   Point(2) = {lx, 0, 0, m1};
Point(3) = {lx, ly, 0, m1}; Point(4) = {0, ly, 0, m1};
// Create points for the interface of the first hole
Point(5) = {cntx, cnty1, 0, m2};        Point(6) = {cntx+cntR, cnty1, 0, m2};
Point(7) = {cntx, cnty1+cntR, 0, m2};   Point(8) = {cntx-cntR, cnty1, 0, m2};
Point(9) = {cntx, cnty1-cntR, 0, m2};
// Create points for the interface of the second hole
Point(10) = {cntx, cnty2, 0, m2};        Point(11) = {cntx+cntR, cnty2, 0, m2};
Point(12) = {cntx, cnty2+cntR, 0, m2};   Point(13) = {cntx-cntR, cnty2, 0, m2};
Point(14) = {cntx, cnty2-cntR, 0, m2};
// Create points for the interface of the third hole
Point(15) = {cntx, cnty3, 0, m2};        Point(16) = {cntx+cntR, cnty3, 0, m2};
Point(17) = {cntx, cnty3+cntR, 0, m2};   Point(18) = {cntx-cntR, cnty3, 0, m2};
Point(19) = {cntx, cnty3-cntR, 0, m2};

//***** Create lines **/
// Creat lines for the out boundary
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 4};  Line(4) = {4, 1};
// Create lines for the interface of first hole
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};
// Create lines for the interface of the second hole
Circle(9) = {11, 10, 12};
Circle(10) = {12, 10, 13};
Circle(11) = {13, 10, 14};
Circle(12) = {14, 10, 11};
// Create lines for the interface of the third hole
Circle(13) = {16, 15, 17};
Circle(14) = {17, 15, 18};
Circle(15) = {18, 15, 19};
Circle(16) = {19, 15, 16};

//*** Create loops from lines **/
Curve Loop(17) = {1, 2, 3, 4};   // Outer boundary
Curve Loop(18) = {6, 7, 8, 5};   // Interface of the first hole
Curve Loop(19) = {10, 11, 12, 9};   // Interface of the second hole
Curve Loop(20) = {14, 15, 16, 13};   // Interface of the third hole

//*** Create blocks from lines **/
Plane Surface(1) = {17, 18, 19, 20}; // creat block SE
Plane Surface(2) = {18}; // creat the first Li block
Plane Surface(3) = {19}; // creat the second Li block

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("right") = {2};
//+
Physical Curve("top") = {3};
//+
Physical Curve("left") = {4};
//+
Physical Curve("pore") = {5:16};
//+
//Physical Curve("interface2") = {9,10,11,12};
//+
//Physical Curve("interface3") = {13,14,15,16};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockLi") = {2,3};
//+
//Physical Surface("blockLi2") = {3};
