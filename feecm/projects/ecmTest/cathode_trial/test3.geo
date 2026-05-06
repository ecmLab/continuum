// Gmsh project created on Mon Nov 16 11:31:14 2020
//+
Point(1) = {0, 0, 0, 20.0};
//+
Point(2) = {4000, 0, 0, 20.0};
//+
Point(3) = {0, -26, 0, 20};
//+
Point(4) = {4000, -26, 0, 20};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Point(5) = {0, 0, 0, 20.0};
//+
Point(6) = {4000, 0, 0, 10.0};
//+
Point(7) = {13000, 0, 0, 50};
//+
Point(8) = {13000, 200, 0, 100};
//+
Point(9) = {0, 200, 0, 20};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 5};
//+
Curve Loop(2) = {5, 6, 7, 8, 9};
//+
Plane Surface(2) = {2};
//+
Point(10) = {0, 200, 0, 20};
//+
Point(11) = {11000, 200, 0, 20};
//+
Point(12) = {11000, 240, 0, 20};
//+
Point(13) = {0, 240, 0, 20};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 13};
//+
Line(13) = {13, 10};
//+
Curve Loop(3) = {10, 11, 12, 13};
//+
Physical Curve("Anode_bottom") =  {10};
//+
Physical Curve("Anode_right") =  {11};
//+
Physical Curve("Anode_top") =  {12};
//+
Physical Curve("Anode_left") =  {13};
//+
Physical Surface("Cathode") = {1};
//+
Physical Surface("Electrolyte") = {2};
//+
Physical Curve("Cathode_bottom") = {3};
//+
Physical Curve("Cathode_right") = {2};
//+
Physical Curve("Cathode_left") = {4};
//+
Physical Curve("Electrolyte_left") = {9};
//+
Physical Curve("Electrolyte_top") = {8};
//+
Physical Curve("Electrolyte_right") = {7};
//+
Physical Curve("Cathode_top") = {1};
//+
Physical Curve("Electrolyte_bottom") = {5, 6};
//+
//+
Plane Surface(3) = {3};
//+
Physical Surface("Anode") = {3};

Transfinite Line
Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
