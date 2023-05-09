/// Model 1: Build a model with a Gr and one Li metal

// Parameters 
xGr   = 26.55;  // thickness of the Gr, in unit um
yGr   = 10.0;  //weight of the Gr, in unit um
xLi   = 0.5;    // depth of the extruded Li, in unit um
yLi   = 0.5;    // width of the extruded Li, in unit um

m1    = 0.1; // mesh characteristic length for Gr-Li interface points
m2    = 0.5; // mesh characteristic length for Gr-SE interface points

//**** Create points **/
Point(1) = {0, 0, 0, m1};     Point(2) = {xGr, 0, 0, m2};
Point(3) = {xGr, yGr, 0, m1}; Point(4) = {xLi, yGr, 0, m1};
Point(5) = {0,  yGr, 0, m1}; 
Point(6) = {xGr, yGr+yLi, 0, m1}; Point(7) = {0, yGr+yLi, 0, m1};

//***** Create lines **/
//AgC
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 4};  Line(4) = {4, 5};
Line(5) = {5, 1};
//Li
Line(6) = {3, 6};  Line(7) = {6, 7};
Line(8) = {7, 5};

//*** Create loops**/
Curve Loop(11)  = {1, 2, 3, 4, 5};   // Gr boundary
//*** Create loops for Li **/
Curve Loop(12)  = {-4, -3, 6, 7, 8};   // Li boundary

//*** Create blocks **/
Plane Surface(1) = {11}; // creat block Gr
Plane Surface(2) = {12}; // creat block Li

//*** Mesh control **/
//Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
//Recombine Surface{1,2};

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("Gr_btm") = {1};
//+
Physical Curve("Gr_rgt") = {2};
//+
Physical Curve("Gr_top") = {3};
//+
Physical Curve("Li_top") = {4};
//+
Physical Curve("Li_lft") = {5};
//+
Physical Curve("Pl_lft") = {8};
//+
Physical Surface("blockGr") = {1};
//+
Physical Surface("blockLi") = {2};
