lx   = 2;       // width of the model, in unit um
ly   = 2;       // height of the model, in unit um
dw   = 0.2;       // width of the defect located at the top middle, in unit um
dh   = 0.1;       // length of the defect, in unit um
dInt = 0.05;      // Initial thickenss of the interlayer, in unit um
nSE  = 21;      // Number of discretization points of the cosine shape defect for SE
nLi  = 31;      // Number of discretization points of the cosine shape defect for Li
nInt = 31;      // Number of discretization points of the cosine shape defect for interlyer

nSE_mesh = 51; 
nLi_mesh = 81;

// Create SE geometry


m0SE  = 0.05;    // mesh characteristic length for bottom points of SE
m1SE  = 0.01;   // mesh characteristic length for interface points of SE
m0Li  = 0.04;    // mesh characteristic length for top points of Li
m1Li  = 0.002;   // mesh characteristic length for bottom points of Li
m0Int = 0.002;   // mesh characteristic length for top points of interlayer
m1Int = 0.002;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
For i In {0 : nSE+1}
    x = dw - dw*i/(nSE + 1);
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
Transfinite Line {2, -4}  = 10 Using Progression 1.1;
Curve Loop(5) = {1:4};
Plane Surface(1) = {5};
Transfinite Surface {1};
Recombine Surface {1};

// *** InterLayer *** // 
For i In {0 : nLi+1}
    x = dw - dw*i/(nLi + 1);
    np= newp;
   intpListTop[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)) + dInt, 0, m1SE};
// Point(np)= {x, dInt, 0, m1SE};
EndFor
For i In {0 : nLi+1}
    x = dw - dw*i/(nLi + 1);
    np= newp;
   intpListBot[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
//    Point(np)= {x, 0, 0, m1SE};
EndFor
Spline(6) = intpListTop[];
Line(7) = {intpListTop[nLi+1], intpListBot[nLi+1]};
Spline(8) = intpListBot[];
Line(9) = {intpListBot[0], intpListTop[0]};
Transfinite Line {6, 8} = nLi_mesh Using Progression 1;
Transfinite Line {7, 9}  = 9 Using Progression 1;
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
Transfinite Line {11, -13} = 10 Using Progression 1.1;
Transfinite Surface {3};
Recombine Surface {3};
//**** */
// ** Rhs of SE ** //
//+
Point(94) = {lx/2, 0.0, 0, m0SE};
//+
Point(95) = {lx/2, -dw, 0, m0SE};
//+
Point(96) = {lx/2, -ly/2, 0, m0SE};
//+
Point(97) = {dw, -ly/2, 0, m0SE};
//+
Point(98) = {0, -ly/2, 0, m0SE};
//+
Line(14) = {1, 94};
//+
Line(15) = {94, 95};
//+
Line(16) = {95, 25};
//+
Line(17) = {95, 96};
//+
Line(18) = {96, 97};
//+
Line(19) = {97, 25};
//+
Line(20) = {97, 98};
//+
Line(21) = {98, 24};
//+
Curve Loop(15) = {4, 14, 15, 16};
//+
Plane Surface(4) = {15};
//+
Curve Loop(16) = {16, -19, -18, -17};
//+
Plane Surface(5) = {16};
//+
Curve Loop(17) = {3, -19, 20, 21};
//+
Plane Surface(6) = {17};
//+
Physical Surface("blockCeramic") = {1, 4, 5, 6};
Transfinite Line {-14, 16, 18} = 51 Using Progression 0.91;
Transfinite Line {17, -19, -21} = 21 Using Progression 1.2;
Transfinite Line{15} = 10 Using Progression 1.1;
Transfinite Line{20} = nSE_mesh;
Transfinite Surface{1:6};
Recombine Surface{1:6};
/// *** Rhs of interLayer *** // 
//+
Point(99) = {lx/2, 0, 0, m0SE};
//+
Point(100) = {lx/2, dInt, 0, m0SE};
// + 
Line(22) = {intpListBot[0], 99};
//+
Line(23) = {intpListTop[0], 100};
//+
Line(24) = {99, 100};
//+
Curve Loop(18) = {24, -23, -9, 22} ;
Plane Surface(7) = {18};//+
Physical Surface("interLayer") = {7, 2};
Transfinite Line {24} = 9 Using Progression 1;
Transfinite Line {22, 23} = 50 Using Progression 1.1;
Transfinite Surface{7,8};
Recombine Surface{7,8};
// *** Metal Block Surface **//
//+
Point(101) = {lx/2, dw, 0, m0SE};
//+
Point(102) = {lx/2, ly/4, 0, m0SE};
//+
Point(103) = {dw, ly/4, 0, m0SE};
//+
Point(104) = {0, ly/4, 0, m0SE};
//+
Line(25) = {100, 101};
//+
Line(26) = {93, 101};
//+
Line(27) = {93, 103};
//+
Line(28) = {103, 102};
//+
Line(29) = {101, 102};
//+
Line(30) = {104, 103};
//+
Line(31) = {92, 104};
//+
Curve Loop(19) = {25, -26, 13, 23};
//+
Plane Surface(8) = {19};
//+
Curve Loop(20) = {29, -28, -27, 26};
//+
Plane Surface(9) = {20};
//+
Curve Loop(21) = {27, -30, -31, 12};
//+
Plane Surface(10) = {21};
//+
Physical Surface("blockMetal") = {8, 9, 10, 3};
//+
Transfinite Line{26,28} = 50 Using Progression 1.1;
//+
Transfinite Line{25} = 10 Using Progression 1.1;
//+
Transfinite Line{31,27, 29} = 20 Using Progression 1.1;
//+
Transfinite Line{30} = nLi_mesh Using Progression 1;
Transfinite Surface{8:10};
Recombine Surface{8:10};
// Compound Curve {20, 18};

// //+
Physical Curve("blockCeramic_left") = {21, 2};
// //+
Physical Curve("blockCeramic_right") = {15, 17};
// //+
Physical Curve("blockMetal_right") = {24, 25, 29};
// //+
Physical Curve("blockMetal_left") = {7, 11, 31};
// //+
Physical Curve("blockMetal_top") = {30, 28};
// //+
Physical Curve("blockMetal_bottom") = {8, 22};
// //+
Physical Curve("blockCeramic_bottom") = {18, 20};
// //+
Physical Curve("blockCeramic_top") = {1, 14};



// //+
// Point(105) = {0, -ly/2 - dInt, 0, 1.0};
// //+
// Point(106) = {lx/2, -ly/2 - dInt, 0, 1.0};
// //+
// Line(32) = {98, 105};
// //+
// Line(33) = {105, 106};
// //+
// Line(34) = {106, 96};
// //+
// Curve Loop(22) = {33, 34, 18, 20, 32};
// //+
// Plane Surface(11) = {22};
// //+
// Physical Surface("blockCeramic") = {11, 5, 4, 1, 6};
// //+
// Physical Curve("blockCeramic_left") = {2, 21, 32};
// //+
// Physical Curve("blockCeramic_right") = {15, 17, 34};
// //+
// Physical Curve("blockCeramic_bottom") = {33};
// //+
// Transfinite Line {33} = 201 Using Progression 1;
// Transfinite Line {32,34} = 6 Using Progression 1;
// // Transfinite Surface{11};
// Recombine Surface{1:11};
//+

