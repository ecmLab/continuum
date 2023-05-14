//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 0.05, 0};
//+
Rectangle(2) = {0, 0.05, 0, 1, 1, 0};
Coherence;
//+
Rectangle(3) = {0, 0, 0, 1, -1, 0};
//+
Physical Surface("interLayer") = {1};
//+
Physical Surface("blockMetal") = {2};
//+
Physical Surface("blockCeramic") = {3};
//+
Physical Curve("blockCeramic_bottom") = {10};
//+
Physical Curve("blockCeramic_right") = {9};
//+
Physical Curve("blockCeramic_left") = {11};
//+
Physical Curve("blockCeramic_top") = {8};
//+
Physical Curve("blockMetal_bottom") = {1};
//+
Physical Curve("blockMetal_left") = {7, 4};
//+
Physical Curve("blockMetal_right") = {5, 2};
//+
Physical Curve("blockMetal_top") = {6};
//+
Transfinite Curve {2, 4} = 3 Using Progression 1;
//+
Transfinite Curve {7, 5} = 5 Using Progression 1;
//+
Transfinite Curve {11, 9} = 5 Using Progression 1;
//+
Transfinite Curve {6, 3, 1} = 10 Using Progression 1;
//+
Transfinite Curve {10, 8} = 6 Using Progression 1;
//+
Transfinite Surface {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Recombine Surface {3, 1, 2};
