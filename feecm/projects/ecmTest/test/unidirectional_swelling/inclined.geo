// Gmsh project created on Fri Aug 14 11:19:34 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1.0, 0.58, 0, 1.0};
//+
Point(3) = {1.0, 0.68, 0, 1.0};
//+
Point(4) = {0, 0.1, 0, 1.0};
//+
Point(5) = {0, 1.0, 0, 1.0};
//+
Point(6) = {1.0, 1.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 3};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 3};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {5} = 4 Using Progression 1;
//+
Transfinite Curve {4, 2} = 2 Using Progression 1;
//+
Transfinite Curve {1, 3, 6} = 5 Using Progression 1;
//+
Physical Surface("interLayer") = {1};
//+
Physical Surface("metalLayer") = {2};
//+
Physical Curve("top") = {6};
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("left") = {5, 4};
//+
Physical Curve("right") = {7, 2};
//+
Recombine Surface {2};
//+
Recombine Surface {1};
