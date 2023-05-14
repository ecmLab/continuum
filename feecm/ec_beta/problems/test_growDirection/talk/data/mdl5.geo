/// Model 1: Build a model with a SE and one Li metal

// Parameters 
xSE   = 400.0;  // length of the SE, in unit um
ySE   = 400.0;  //height of the SE, in unit um
xLi   = 100;    // length of the Li, in unit um
yLi   = 1;     //height of the Li, in unit um
rp    = 0.1;  // the round radius at corner point, in unit um

m1    = 10; // mesh characteristic length for boundary points
m2    = 0.1; // mesh characteristic length for interface points
m3    = 0.01; // mesh characteristic length for inner points

//**** Create points **/
Point(1) = {0, 0, 0, m1};     Point(2) = {xSE, 0, 0, m1};
Point(3) = {xSE, ySE, 0, m1}; Point(4) = {xLi+rp, ySE, 0, m3};
Point(5) = {rp,  ySE, 0, m3}; Point(6) = {0, ySE-rp, 0, m3}; 
Point(7) = {-yLi, ySE-rp, 0, m2}; Point(8) = {rp, ySE+yLi, 0, m2}; 
Point(9) = {xLi+rp, ySE+yLi, 0, m2}; Point(10) = {rp, ySE-rp, 0, m3}; 
Point(11) = {xLi+rp, ySE+yLi/2, 0, m3}; 

//***** Create lines **/
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 4};  Line(4) = {4, 5};
Circle(5) = {5, 10, 6};  Line(6) = {6, 1};  
Line(7) = {6, 7};  Circle(8) = {7, 10, 8};
Line(9) = {8, 9};  Circle(10) = {4, 11, 9};

//*** Create loops from lines **/
Curve Loop(11)  = {1, 2, 3, 4, 5, 6};   // SE boundary
Curve Loop(12)  = {4, 5, 7, 8, 9, -10};   // Li boundary

//*** Create blocks from lines **/
Plane Surface(1) = {11}; // creat block SE
Plane Surface(2) = {12}; // creat block Li

//*** Mesh control **/
//Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
Recombine Surface{1,2};

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("SE_bottom") = {1};
//+
Physical Curve("SE_right") = {2};
//+
Physical Curve("SE_top") = {3};
//+
Physical Curve("SE_left") = {6};
//+
Physical Curve("Li_right") = {10};
//+
Physical Curve("Li_top") = {8,9};
//+
Physical Curve("Li_left") = {7};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockLi") = {2};
