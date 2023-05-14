// Gmsh project created on Mon Nov 16 11:31:14 2020
//+
Point(1) = {0, 0, 0, 40.0};
//+
Point(2) = {4000, 0, 0, 40};
//+
Point(3) = {0, -26, 0, 40};
//+
Point(4) = {4000, -26, 0, 40};
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
Point(5) = {0, 0, 0, 40.0};
//+
Point(6) = {4000, 0, 0, 20};
//+
Point(7) = {13000, 0, 0, 100};
//+
Point(8) = {13000, 200, 0, 100};
//+
Point(9) = {0, 200, 0, 100};
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