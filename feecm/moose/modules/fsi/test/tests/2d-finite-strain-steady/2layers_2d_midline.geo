SetFactory ("OpenCASCADE");
r1 = 0.4;
r2 = 0.8;
height = 10;
lc = 1;

Point(1) = {0, 0, 0, lc};
Point(2) = {0, height, 0, lc};
Point(3) = {r1, 0, 0, lc};
Point(4) = {r1, height, 0, lc};
Point(5) = {-r1, 0, 0, lc};
Point(6) = {-r1, height, 0, lc};
Point(7) = {r2, 0, 0, lc};
Point(8) = {r2, height, 0, lc};
Point(9) = {-r2, 0, 0, lc};
Point(10) = {-r2, height, 0, lc};

Line(1) = {1, 2};
Transfinite Line{1} = 6;
Line(2) = {2, 4};
Transfinite Line{2} = 4;
Line(3) = {4, 3};
Transfinite Line{3} = 6;
Line(4) = {3, 1};
Transfinite Line{4} = 4;
Line(5) = {2, 6};
Transfinite Line{5} = 4;
Line(6) = {6, 5};
Transfinite Line{6} = 6;
Line(7) = {5, 1};
Transfinite Line{7} = 4;
Line(8) = {4, 8};
Transfinite Line{8} = 2;
Line(9) = {8, 7};
Transfinite Line{9} = 6;
Line(10) = {7, 3};
Transfinite Line{10} = 2;
Line(11) = {6, 10};
Transfinite Line{11} = 2;
Line(12) = {10, 9};
Transfinite Line{12} = 6;
Line(13) = {9, 5};
Transfinite Line{13} = 2;

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {1, 5, 6, 7};
Curve Loop(3) = {-3, 8, 9, 10};
Curve Loop(4) = {-6, 11, 12, 13};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Surface{4};
Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
Recombine Surface{4};

Physical Point("pin1") = {1};
Physical Point("pin2") = {2};
Physical Curve("fluid_bottom") = {10, 13};
Physical Curve("solid_bottom") = {4, 7};
Physical Curve("fluid_top") = {8, 11};
Physical Curve("solid_top") = {2, 5};
Physical Curve("fluid_wall") = {9, 12};
Physical Curve("solid_wall") = {3, 6};
Physical Surface("fluid") = {3, 4};
Physical Surface("solid") = {1, 2};