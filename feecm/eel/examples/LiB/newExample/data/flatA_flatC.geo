//Geometry of the battery model, all in unit mm
W = 0.01;       // Width of the battery
La = 0.02;      // thickness of the anode
Le = 0.04;      // thickness of the electrolyte
Lc = 0.02;      // thickness of the cathode
nSE  = 101;     // discretization points of the cosine shape defect

e1 = 0.001;     // coarse mesh control
e2 = 0.001;     // fine mesh control

//**** Create conner points of each component**/
Point(1) = {0, 0, 0, e1};
Point(2) = {0, W, 0, e1};
Point(3) = {La, W, 0, e2};
Point(4) = {La+Le, W, 0, e2};
Point(5) = {La+Le+Lc, W, 0, e1};
Point(6) = {La+Le+Lc, 0, 0, e1};
Point(7) = {La+Le, 0, 0, e2};
Point(8) = {La, 0, 0, e2};
//**** Create lines of each component**/
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {3, 8};
Line(10) = {4, 7};

Line Loop(1) = {1, 2, 9, 8};
Line Loop(2) = {-9, 3, 10, 7};
Line Loop(3) = {-10, 4, 5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

//Control of mesh quality
Mesh.CharacteristicLengthMin = 0.0001;
Mesh.CharacteristicLengthMax = 0.003;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.Smoothing = 10;
Mesh.Optimize = 1;
Recombine Surface{1,2,3};

Physical Line("left") = {1};
Physical Line("right") = {5};
Physical Line("top") = {2, 3, 4};
Physical Line("bottom") = {6, 7, 8};
Physical Surface("anode") = {1};
Physical Surface("elyte") = {2};
Physical Surface("cathode") = {3};
