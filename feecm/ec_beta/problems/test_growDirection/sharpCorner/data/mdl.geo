/// Model 1: Build a model with a SE and one Li metal

// Parameters 
xSE   = 0.5;  // length of the SE, in unit um
ySE   = 0.5;  //height of the SE, in unit um
xLi   = 0.02;  // length of the Li, in unit um
yLi   = 0.001;  //height of the Li, in unit um

m1    = 0.02; // mesh characteristic length for boundary points
m2    = 0.001; // mesh characteristic length for interface points
m3    = 0.00002; // mesh characteristic length for inner points

//**** Create points **/
Point(1) = {0, 0, 0, m1};     Point(2) = {xSE, 0, 0, m1};
Point(3) = {xSE, ySE, 0, m1}; Point(4) = {0, ySE, 0, m3};
Point(5) = {xLi, ySE, 0, m3}; Point(6) = {xLi, ySE+yLi, 0, m3}; 
Point(7) = {0, ySE+yLi, 0, m3};

//***** Create lines **/
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 5};  Line(4) = {5, 4};
Line(5) = {4, 1};  Line(6) = {5, 6};  
Line(7) = {6, 7};  Line(8) = {7, 4};

//*** Create loops from lines **/
Curve Loop(10)  = {1, 2, 3, 4, 5};   // SE boundary
Curve Loop(11)  = {-4,6, 7, 8};   // Li boundary

//*** Create blocks from lines **/
Plane Surface(1) = {10}; // creat block SE
Plane Surface(2) = {11}; // creat block Li

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("SE_bottom") = {1};
//+
Physical Curve("SE_right") = {2};
//+
Physical Curve("SE_left") = {5};
//+
Physical Curve("Li_right") = {6};
//+
Physical Curve("Li_top") = {7};
//+
Physical Curve("Li_left") = {8};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockLi") = {2};
