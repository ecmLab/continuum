lx = 0.5;
ly_si = 0.1;
ly_se = 0.5;
nx_si = 50;
nx_se = 40;
ny_si = 20;
ny_se = 20;

// Gmsh project created on Mon Nov 16 11:31:14 2020
//+
Point(1) = {0, 0, 0, 20.0};
//+
Point(2) = {lx, 0, 0, 20.0};
//+
Point(3) = {0, -ly_se, 0, 20};
//+
Point(4) = {lx, -ly_se, 0, 20};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Point(5) = {0, 0, 0, 40.0};
//+
Point(6) = {lx, 0, 0, 50};
//+
Point(7) = {lx, ly_si, 0, 40};
//+
Point(8) = {0, ly_si, 0, 40};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Electrolyte") = {1};
//+
Physical Surface("Anode") = {2};
//+
Physical Curve("Electrolyte_top") = {1};
//+
Physical Curve("Electrolyte_right") = {2};
//+
Physical Curve("Electrolyte_bottom") = {3};
//+
Physical Curve("Electrolyte_left") = {4};
Transfinite Curve {2, 4} = ny_se Using Progression 1;
Transfinite Curve {1, 3} = nx_se Using Progression 1;
//+
Physical Curve("Anode_bottom") = {5};
//+
Physical Curve("Anode_right") = {6};
//+
Physical Curve("Anode_top") = {7};
//+
Physical Curve("Anode_left") = {8};
//+
Transfinite Curve {6, 8} = ny_si Using Progression 1;
Transfinite Curve {5, 7} = nx_si Using Progression 1;
//+
Transfinite Surface {2};
//+
Transfinite Surface {1};
//+
Recombine Surface{1, 2};
