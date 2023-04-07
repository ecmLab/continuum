//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 400, 0.5, 0};
//+
Rectangle(2) = {0, 0.5, 0, 1000, 200, 0};
//+
Physical Surface("electrolyte") = {2};
//+
Physical Surface("cathode") = {1};
//+
Characteristic Length {3} = 0.05;
//+
Transfinite Curve {4, 2} = 5 Using Progression 1;
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {2};
