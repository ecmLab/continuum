lx   = 2;       // width of the model, in unit um
ly   = 2;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.05;       // length of the defect, in unit um
dInt = 0.01;      // Initial thickenss of the interlayer, in unit um
nSE  = 21;      // Number of discretization points of the cosine shape defect for SE
nLi  = 31;      // Number of discretization points of the cosine shape defect for Li
nInt = 31;      // Number of discretization points of the cosine shape defect for interlyer

nSE_mesh = 41; 
nLi_mesh = 61;

// Create SE geometry


m0SE  = 0.05;    // mesh characteristic length for bottom points of SE
m1SE  = 0.01;   // mesh characteristic length for interface points of SE
m0Li  = 0.04;    // mesh characteristic length for top points of Li
m1Li  = 0.002;   // mesh characteristic length for bottom points of Li
m0Int = 0.002;   // mesh characteristic length for top points of interlayer
m1Int = 0.002;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
For i In {0 : nSE+1}
    x = dw - 2*dw*i/(nSE + 1);
    np= newp;
   SEpList[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
EndFor
Spline(1) = SEpList[];
coord = Point{SEpList[nSE + 1]};
np  = newp;
Point(np) = {coord[0], -dw, 0, m0SE};
Line(2) = {SEpList[nSE + 1], np};
coord = Point{SEpList[0]};
np = newp;
Point(np) = {coord[0], - dw, 0, m0SE};
Line(3) = {np-1, np};
Line(4) ={np, SEpList[0]};
Transfinite Line {1, 3} = nSE_mesh Using Progression 1;
Transfinite Line {2, -4}  = 5 Using Progression 1.1;
Curve Loop(5) = {1:4};
Plane Surface(1) = {5};
Transfinite Surface {1};
Recombine Surface {1};

// *** InterLayer *** // 
For i In {0 : nLi+1}
    x = dw - 2*dw*i/(nLi + 1);
    np= newp;
   intpListTop[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)) + dInt, 0, m1SE};
EndFor
For i In {0 : nLi+1}
    x = dw - 2*dw*i/(nLi + 1);
    np= newp;
   intpListBot[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
EndFor
Spline(6) = intpListTop[];
Line(7) = {intpListTop[nLi+1], intpListBot[nLi+1]};
Spline(8) = intpListBot[];
Line(9) = {intpListBot[0], intpListTop[0]};
Transfinite Line {6, 8} = nLi_mesh Using Progression 1;
Transfinite Line {7, 9}  = 5 Using Progression 1;
Curve Loop(10) = {6, 7, -8, 9};
Plane Surface(2) = {10};
Transfinite Surface {2};
Recombine Surface {2};

// *** Metal **** // 
// Create box for transfinite
np = newp;
coord = Point{intpListTop[nLi + 1]};
Point(np) = { coord[0], dw, 0, m0SE}; 
Line(11) = {intpListTop[nLi + 1], np};
np = newp;
coord = Point{intpListTop[0]};
Point(np) = { coord[0], dw, 0, m0SE}; 
Line(12) = {np-1, np};
Line(13) = {np, intpListTop[0]};
Curve Loop(14) = {6, 11, 12, 13};
Plane Surface(3) = {14};
Transfinite Line {12} = nLi_mesh Using Progression 1;
Transfinite Line {11, -13} = 5 Using Progression 1.1;
Transfinite Surface {3};
Recombine Surface {3};
//+
//+
Point(94) = {-lx/2, 0, 0, m0SE};
//+
Point(95) = {-lx/2, -dw, 0, m0SE};
//+
Point(96) = {-lx/2, -ly/2, 0, m0SE};
//+
Point(97) = {-dw, -ly/2, 0, m0SE};
//+
Point(98) = {dw, -ly/2, 0, m0SE};
//+
Point(99) = {lx/2, -ly/2, 0, m0SE};
//+
Point(100) = {lx/2, -dw, 0, m0SE};
//+
Point(101) = {lx/2, 0, 0, m0SE};
//+
Line(14) = {23, 94};
//+
Line(15) = {94, 95};
//+
Line(16) = {95, 24};
//+
Line(17) = {95, 96};
//+
Line(18) = {96, 97};
//+
Line(19) = {97, 24};
//+
Line(20) = {97, 98};
//+
Line(21) = {98, 25};
//+
Line(22) = {98, 99};
//+
Line(23) = {99, 100};
//+
Line(24) = {100, 25};
//+
Line(25) = {100, 101};
//+
Line(26) = {101, 1};
//+
Point(102) = {-lx/2, dInt, 0, m0SE};
//+
Point(103) = {lx/2, dInt, 0, m0SE};
//+
Line(27) = {58, 102};
//+
Line(28) = {103, 26};
//+
Point(104) = {-lx/2, 0, 0, m0SE};
//+
Point(105) = {lx/2, 0, 0, m0SE};
//+
Line(29) = {intpListBot[nLi+1], 104};
//+
Line(30) = {105, intpListBot[0]};//+
Point(106) = {-lx/2, dw, 0, m0SE};
//+
Point(107) = {lx/2, dw, 0, m0SE};
//+
Line(31) = {92, 106};
//+
Line(32) = {106, 102};
//+
Line(33) = {103, 107};
//+
Line(34) = {107, 93};
//+
Line(35)  = {105,103};
//+
Line(36) = {102, 104};//+
Point(108) = {-lx/2, ly/4, -0, m0SE};
//+
Point(109) = {-dw, ly/4, -0, m0SE};
//+
Point(110) = {dw, ly/4, -0, m0SE};
//+
Point(111) = {lx/2, ly/4, -0, m0SE};
//+
Line(37) = {107, 111};
//+
Line(38) = {111, 110};
//+
Line(39) = {110, 93};
//+
Line(40) = {110, 109};
//+
Line(41) = {109, 92};
//+
Line(42) = {109, 108};
//+
Line(43) = {108, 106};
//+
Curve Loop(15) = {2, -16, -15, -14};
//+
Plane Surface(4) = {15};
//+
Curve Loop(16) = {19, -16, 17, 18};
//+
Plane Surface(5) = {16};
//+
Curve Loop(17) = {3, -21, -20, 19};
//+
Plane Surface(6) = {17};
//+
Curve Loop(18) = {24, -21, 22, 23};
//+
Plane Surface(7) = {18};
//+
Curve Loop(19) = {25, 26, -4, -24};
//+
Plane Surface(8) = {19};
//+
Physical Surface("blockCeramic") = {1, 4, 5, 6, 7, 8};
//+
Curve Loop(20) = {11, 31, 32, -27};
//+
Plane Surface(9) = {20};
//+
Curve Loop(21) = {41, 31, -43, -42};
//+
Plane Surface(10) = {21};
//+
Curve Loop(22) = {12, -39, 40, 41};
//+
Plane Surface(11) = {22};
//+
Curve Loop(23) = {37, 38, 39, -34};
//+
Plane Surface(12) = {23};
//+
Curve Loop(24) = {33, 34, 13, -28};
//+
Plane Surface(13) = {24};
//+
Physical Surface("blockMetal") = {3, 13, 12, 11, 10, 9};
//+
Curve Loop(25) = {7, 29, -36, -27};
//+
Plane Surface(14) = {25};
//+
Curve Loop(26) = {28, -9, -30, 35};
//+
Plane Surface(15) = {26};
//+
Physical Surface("interLayer") = {15, 2, 14};
//+
Physical Curve("blockCeramic_left") = {17, 15};
//+
Physical Curve("blockCeramic_right") = {25, 23};
//+
Physical Curve("blockCeramic_bottom") = {18, 20, 22};
//+
Physical Curve("blockCeramic_top") = {14, 1, 26};
//+
Physical Curve("blockMetal_top") = {42, 40, 38};
//+
Physical Curve("blockMetal_left") = {36, 32, 43};
//+
Physical Curve("blockMetal_right") = {35, 33, 37};
//+
Physical Curve("blockMetal_bottom") = {29, 8, 30};
//+

// Transfinite Line{20} = nSE Using Progression 1;
Transfinite Line {36, 35} = 5 Using Progression 1;

Transfinite Line {40} = nLi_mesh Using Progression 1;

Transfinite Line {15, -25} = 5 Using Progression 1.1;
Transfinite Line {20} = nSE_mesh Using Progression 1;

Transfinite Line {18, 16, 14} = 50 Using Progression 1;
Transfinite Line {22, 24, 26} = 50 Using Progression 1;

Transfinite Line {30, 28, 34, 38} = 50 Using Progression 0.91;
Transfinite Line {29, 27, 31, 42} = 50 Using Progression 1.1;
Transfinite Line {43, 41, 39, 37} = 15 Using Progression 1;

Transfinite Line {32, 33} = 5 Using Progression 1.1;

Transfinite Line {17, -19} = 15 Using Progression 1.1;
Transfinite Line {21, 23} = 15 Using Progression 0.91;
Transfinite Surface{1:18};

Recombine Surface{1:18};
