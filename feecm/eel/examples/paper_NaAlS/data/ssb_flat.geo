Mesh.CharacteristicLengthMin = 0.001;
Mesh.CharacteristicLengthMax = 0.003;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.Smoothing = 10;
Mesh.Optimize = 1;
W = 0.01;
Lc = 0.02;
Le = 0.02;
La = 0.01;
R  = 0.0025;
e1 = 0.001;
e2 = 0.0001;

Point(1) = {0, 0, 0, e1};
Point(2) = {0, W, 0, e1};
Point(3) = {Lc, W, 0, e1};
Point(4) = {Lc+Le, W, 0, e2};
Point(5) = {Lc+Le+La, W, 0, e1};
Point(6) = {Lc+Le+La, 0, 0, e1};
Point(7) = {Lc+Le, 0, 0, e2};
Point(8) = {Lc, 0, 0, e1};

Point(9) = {Lc/2, W/2, 0, e2};
Point(10) = {Lc/2+R, W/2, 0, e2};
Point(11) = {Lc/2, W/2+R, 0, e2};
Point(12) = {Lc/2-R, W/2, 0, e2};
Point(13) = {Lc/2, W/2-R, 0, e2};

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

Circle(11) = {10, 9, 11};
Circle(12) = {11, 9, 12};
Circle(13) = {12, 9, 13};
Circle(14) = {13, 9, 10};

Line Loop(1) = {1, 2, 9, 8};
Line Loop(2) = {-9, 3, 10, 7};
Line Loop(3) = {-10, 4, 5, 6};
Line Loop(4) = {11, 12, 13, 14};

Plane Surface(1) = {4};
Plane Surface(2) = {1, 4};
Plane Surface(3) = {2};
Plane Surface(4) = {3};

Physical Line("left") = {1};
Physical Line("right") = {5};
Physical Line("top") = {2, 3, 4};
Physical Line("bottom") = {6, 7, 8};
Physical Surface("cp") = {1};
Physical Surface("cm") = {2};
Physical Surface("e") = {3};
Physical Surface("a") = {4};
