//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, -10, 0, 1, 10, 0};
//+
Rectangle(2) = {0.5, 0, 0, 0.5, 1, 0};
//+
Rectangle(3) = {0.0, 0, 0, 0.5, 0.1, 0};
//+
Rectangle(4) = {0.0, 0.1, 0, 0.5, 0.9, 0};
//+
Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockPore") = {2};
//+
Physical Surface("blockMetal") = {4};
//+
Physical Surface("interLayer") = {3};
//+
Physical Curve("blockCeramic_bottom") = {1};
//+
Physical Curve("blockCeramic_top") = {3};
//+
Physical Curve("blockCeramic_left") = {4};
//+
Physical Curve("blockCeramic_right") = {2};
//+
Physical Curve("blockMetal_top") = {15};
//+
Physical Curve("blockMetal_left") = {16};
//+
Physical Curve("blockMetal_right") = {14};
//+
Physical Curve("blockMetal_bottom") = {13};
//+
Physical Curve("blockPore_left") = {8};
//+
Physical Curve("blockPore_right") = {6};
//+
Physical Curve("blockPore_top") = {5};
//+
Physical Curve("blockPore_bottom") = {7};
//+
Physical Curve("interLayer_left") = {12};
//+
Physical Curve("interLayer_right") = {10};
//+
Physical Curve("interLayer_top") = {11};
//+
Physical Curve("interLayer_bottom") = {9};


//+
Transfinite Curve {15, 13, 11,9, 7, 5, 3, 1} = 15 Using Progression 1;
//+
Transfinite Curve {12, 10} = 20 Using Progression 1;
//+
Transfinite Curve {16, 14, 8,6, 2, 4} = 10 Using Progression 1;


//+
Transfinite Surface {1};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {3};
//+
Recombine Surface {1, 2, 4, 3};
//+
Recombine Surface {1};
