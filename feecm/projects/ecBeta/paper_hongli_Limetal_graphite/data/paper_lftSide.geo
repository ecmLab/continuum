/// Model 1: Build a model with a Gr and one Li metal

// Parameters 
tGr   = 26.55;  // thickness of the Gr, in unit um
wGr   = 10.0;  //weight of the Gr, in unit um
dLi   = 24;    // depth of the extruded Li, in unit um

m1    = 0.1; // mesh characteristic length for Gr-Li interface points
m2    = 0.5; // mesh characteristic length for Gr-SE interface points

//**** Create points **/
Point(1) = {0, 0, 0, m2};         Point(2) = {wGr, 0, 0, m2};
Point(3) = {wGr, tGr-dLi, 0, m1}; Point(4) = {wGr, tGr, 0, m1};
Point(5) = {0,  tGr, 0, m1}; 

//***** Create lines **/
Line(1) = {1, 2};  Line(2) = {2, 3};  
Line(3) = {3, 4};  Line(4) = {4, 5};
Line(5) = {5, 1};  

//*** Create loops for Gr **/
Curve Loop(11)  = {1, 2, 3, 4, 5};   // Gr boundary

//*** Create blocks for Gr **/
Plane Surface(1) = {11}; // creat block Gr

//*** Mesh control **/
//Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
//Recombine Surface{1,2};

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("Gr_SE") = {1};
//+
Physical Curve("Gr_rgt") = {2};
//+
Physical Curve("Gr_lft") = {5};
//+
Physical Curve("Gr_Li") = {3,4};
//+
Physical Surface("blockGr") = {1};
