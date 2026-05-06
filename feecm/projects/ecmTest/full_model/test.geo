//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 1e-2, 0, 1, 1, 0};
//+
Rectangle(2) = {0, 0, 0, 1, -1, 0};
//+
Physical Surface("blockMetal") = {1};
//+
Physical Surface("blockCeramic") = {2};
//+
Physical Curve("blockMetal_left") = {4};
//+
Physical Curve("blockMetal_top") = {3};
//+
Physical Curve("blockMetal_right") = {2};
//+
Physical Curve("blockMetal_bottom") = {1};
//+
Physical Curve("blockCeramic_bottom") = {7};
//+
Physical Curve("blockCeramic_right") = {6};
//+
Physical Curve("blockCeramic_left") = {8};
//+
Physical Curve("blockCeramic_top") = {5};
//+
Transfinite Curve{1,2, 3, 4} = 5 Using Progression 1;
Transfinite Curve{5,6, 7, 8} = 4 Using Progression 1;
Transfinite Surface{1}  = {1, 2, 3, 4};
Transfinite Surface{2}  = {5, 6, 7, 8};
Recombine Surface{1};
Recombine Surface{2};