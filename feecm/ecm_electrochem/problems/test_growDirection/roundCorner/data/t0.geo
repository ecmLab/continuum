/// Build a model with a SE and one Li metal

// Parameters 
tSE   = 20.0;  // thickness of the SE, in unit um
wSE   = 8.0;  // width of the SE, in unit um
tLi   = 0.1;    // initial thickness of the Li, in unit um
//hLi   = 0.5;  // initial height of the Li, in unit um
rp    = 0.05;  // the round radius at corner point, in unit um

m1    = 0.5; // mesh characteristic length for SE boundary points
m2    = 0.05; // mesh characteristic length for Li metal points
m3    = 0.005; // mesh characteristic length for interface points

//**** Create points **/
Point(1) = {0, 0, 0, m1};        Point(2) = {wSE, 0, 0, m1};
Point(3) = {wSE, tSE-rp, 0, m3}; Point(4) = {wSE-rp, tSE-rp, 0, m3};
Point(5) = {wSE-rp, tSE, 0, m3}; Point(6) = {0, tSE, 0, m3}; 
Point(7) = {0, tSE+tLi, 0, m2};  Point(8) = {wSE-rp, tSE+tLi, 0, m2}; 
Point(9) = {wSE+tLi, tSE-rp, 0, m2};

//***** Create lines **/
Line(1)   = {1, 2};     Line(2)   = {2, 3};  
Circle(3) = {3, 4, 5};  Line(4)   = {5, 6};
Line(5)   = {6, 1};     Line(6)   = {6, 7};     
Line(7)   = {7, 8};     Circle(8) = {8, 4, 9};  
Line(9)   = {9, 3};

//*** Create loops from lines **/
Curve Loop(10)  = {1, 2, 3, 4, 5};   // SE boundary
Curve Loop(11)  = {3, 4, 6, 7, 8, 9};   // Li boundary

//*** Create blocks from lines **/
Plane Surface(1) = {10}; // creat block SE
Plane Surface(2) = {11}; // creat block Li

//*** Mesh control **/
//Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
Recombine Surface{1,2};

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("SE_bottom") = {1};
//+
Physical Curve("SE_right")  = {2};
//+
Physical Curve("SE_top")    = {3, 4};
//+
Physical Curve("SE_left")   = {5};
//+
Physical Curve("Li_right")  = {10};
//+
Physical Curve("Li_top")    = {7,8};
//+
Physical Curve("Li_left")   = {6};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockLi") = {2};
