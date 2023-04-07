/// Model 1: Build a model with three void vertically distributed

// Parameters 
lx   = 5;  // length of the model, in unit um
ly   = 5;  //height of the model, in unit um
cntR = 1.25; // Radius of the hole
cntx = lx/2; // x coordinate of the hole center
cnty = ly/2; // y coordinate of the first hole center

m1    = 0.05; // mesh characteristic length for boundary points
m2    = 0.005; // mesh characteristic length for interface points

//**** Create points **/
// Create points for the boundary
Point(1) = {0, 0, 0, m1};   Point(2) = {lx, 0, 0, m1};
Point(3) = {lx, ly, 0, m1}; Point(4) = {0, ly, 0, m1};
// Create points for the interface of the first hole
Point(5) = {cntx, cnty, 0, m2};        Point(6) = {cntx+cntR, cnty, 0, m2};
Point(7) = {cntx, cnty+cntR, 0, m2};   Point(8) = {cntx-cntR, cnty, 0, m2};
Point(9) = {cntx, cnty-cntR, 0, m2};

//***** Create lines **/
// Creat lines for the boundary
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 4};  Line(4) = {4, 1};
// Create lines for the interface of first hole
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

//*** Create loops from lines **/
Curve Loop(17) = {1, 2, 3, 4};   // Outer boundary
Curve Loop(18) = {6, 7, 8, 5};   // Interface of the first hole

//*** Create blocks from lines **/
Plane Surface(1) = {17, 18}; // creat block SE
Plane Surface(2) = {18}; // creat block NMC

