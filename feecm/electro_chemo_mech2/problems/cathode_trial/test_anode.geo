// Gmsh project created on Mon Nov 16 11:31:14 2020
//+
l_cathode = 4000;
l_electrolyte = 13000;
l_anode = 11000;
y_cathode = 26;
y_electrolyte = 200;
y_anode = 40;

Point(1) = {0, 0, 0, 20.0};
//+
Point(2) = {l_cathode, 0, 0, 20.0};
//+
Point(3) = {l_anode, 0, 0, 20.0};
//+
Point(4) = {0, -y_anode, 0, 20};
//+
Point(5) = {l_cathode, -y_anode, 0, 20};
//+
Point(6) = {l_anode, -y_anode, 0, 20.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {1, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 3};
//+
Curve Loop(1) = {-2, -1, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Point(7) = {0, 0, 0, 20.0};
//+
Point(8) = {l_cathode, 0, 0, 10.0};
//+
Point(9) = {l_anode, 0, 0, 50};
//+
Point(10) = {l_electrolyte, 0, 0, 100};
//+
Point(11) = {0, y_electrolyte, 0, 20};
//+
Point(12) = {l_cathode, y_electrolyte, 0, 20};
//+
Point(13) = {l_electrolyte, y_electrolyte, 0, 20};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 13};
//+
Line(11) = {13, 12};
//+
Line(12) = {12, 11};
//+
Line(13) = {11, 7};
//+
Curve Loop(2) = {7,8,9,10,11,12,13};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Anode") = {1};
//+
Physical Surface("Electrolyte") = {2};
//+
Physical Curve("Anode_bottom") = {4,5};
//+
Physical Curve("Anode_right") = {6};
//+
Physical Curve("Anode_left") = {3};
//+
Physical Curve("Electrolyte_left") = {13};
//+
Physical Curve("Electrolyte_top") = {11,12};
//+
Physical Curve("Electrolyte_right") = {10};
//+
Physical Curve("Anode_top") = {1,2};
//+
Physical Curve("Electrolyte_bottom") = {7, 8, 9};
//+
Physical Curve("Electrolyte_current") = {12};
//+
Recombine Surface{1};
Recombine Surface{2}; */
