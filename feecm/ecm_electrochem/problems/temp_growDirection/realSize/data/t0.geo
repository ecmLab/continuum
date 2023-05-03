/// Build a model with a SE and two Li metal as electrodes

// Parameters 
xSE   = 500.0;   // SE thickness, in unit um
ySE   = 200.0;   // SE width, in unit um
xLi   = 100;     // length of the Li, in unit um
yLi   = 5;       // Li initial thickness, in unit um
xCd   = 100;  // Thickness of Cathode, in unit um
rp    = 5;  // the round radius at corner point, in unit um

m1    = 1.0; // mesh characteristic length for interface points
m2    = 2.0; // mesh characteristic length for boundary points of Anode
m3    = 4.0; // mesh characteristic length for boundary points of cathode

//**** Create points **/
Point(1) = {0, 0, 0, m1};     Point(2) = {xSE, 0, 0, m1};
Point(3) = {xSE, ySE, 0, m1}; Point(4) = {xLi+rp, ySE, 0, m1};
Point(5) = {rp,  ySE, 0, m1}; Point(6) = {0, ySE-rp, 0, m1}; 
Point(7) = {-yLi, 0, 0, m2};  Point(8) = {-yLi, ySE+yLi, 0, m2}; 
Point(9) = {xLi+rp, ySE+yLi, 0, m2}; Point(10) = {rp, ySE-rp, 0, m2}; 
Point(11) = {xLi+rp, ySE+yLi/2, 0, m2}; 
Point(12) = {xSE+xCd, 0, 0, m3}; Point(13) = {xSE+xCd, ySE, 0, m3};

//***** Create lines **/
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 4};  Line(4) = {4, 5};
Circle(5) = {5, 10, 6};  
Line(6) = {6, 1};  Line(7) = {1, 7};  
Line(8) = {7, 8};  Line(9) = {8, 9};  
Circle(10) = {4, 11, 9};
Line(11) = {2, 12};  Line(12) = {12, 13};  
Line(13) = {13, 3};  

//*** Create loops from lines **/
Curve Loop(14)  = {1, 2, 3, 4, 5, 6};   // SE boundary
Curve Loop(15)  = {4, 5, 6, 7, 8, 9, -10};   // Anode boundary
Curve Loop(16)  = {11, 12, 13, -2};   // Cathode boundary

//*** Create blocks from lines **/
Plane Surface(1) = {14}; // creat block SE
Plane Surface(2) = {15}; // creat block Anode
Plane Surface(3) = {16}; // creat block Cathode

//*** Mesh control **/
//Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
Recombine Surface{1,2,3};

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("SE_bottom") = {1};
//+
Physical Curve("SE_right") = {2};
//+
Physical Curve("SE_top") = {3,4,5};
//+
//Physical Curve("SE_left") = {6};
//+
Physical Curve("And_right") = {10};
//+
Physical Curve("And_top") = {9};
//+
Physical Curve("And_left") = {8};
//+
Physical Curve("Ctd_bottom") = {11};
//+
Physical Curve("Ctd_right") = {12};
//+
Physical Curve("Ctd_top") = {13};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockAnd") = {2};
//+
Physical Surface("blockCtd") = {3};
