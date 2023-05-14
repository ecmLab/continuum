lx   = 2;       // width of the model, in unit um
ly   = 2;       // height of the model, in unit um
dw   = 0.2;       // width of the defect located at the top middle, in unit um
dh   = 0.4;       // length of the defect, in unit um
dInt = 0.05;      // Initial thickenss of the interlayer, in unit um
nSE  = 21;      // Number of discretization points of the cosine shape defect for SE
nLi  = 31;      // Number of discretization points of the cosine shape defect for Li
nInt = 31;      // Number of discretization points of the cosine shape defect for interlyer

nSE_mesh = 51; 
nLi_mesh = 81;

// Create SE geometry


m0SE  = 0.05;    // mesh characteristic length for bottom points of SE
m1SE  = 0.01;   // mesh characteristic length for interface points of SE
m0Li  = 0.1;    // mesh characteristic length for top points of Li
m1Li  = 0.002;   // mesh characteristic length for bottom points of Li
m0Int = 0.005;   // mesh characteristic length for top points of interlayer
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
Point(np) = {coord[0], -ly/2, 0, m0SE};
Line(2) = {SEpList[nSE + 1], np};
coord = Point{SEpList[0]};
np = newp;
Point(np) = {lx/2, -ly/2, 0, m0SE};
Line(3) = {np-1, np};
np1 = newp;
Point(np1) = {lx/2,0, 0, m0SE};
Line(4) ={np, np1};
Line(5) = {np1, SEpList[0]};
Curve Loop(5) = {1:5};
Plane Surface(1) = {5};
Recombine Surface {1};

// *** InterLayer *** // 
For i In {0 : nLi+1}
    x = dw - dw*i/(nLi + 1);
    np= newp;
   intpListTop[i] = np;
//    Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)) + dInt, 0, m0Int};
   Point(np)= {x, dInt, 0, m0Int};
EndFor
For i In {0 : nLi+1}
    x = dw - dw*i/(nLi + 1);
    np= newp;
   intpListBot[i] = np;
//    Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m0Int};
   Point(np)= {x, 0, 0, m0Int};
EndFor
Spline(6) = intpListTop[];
Line(7) = {intpListTop[nLi+1], intpListBot[nLi+1]};
Spline(8) = intpListBot[];
np = newp;
Point(np) = {lx/2, 0, 0, m0SE};
Line(9) = {np, intpListBot[0]};
np1 = newp;
coord = Point{intpListTop[0]};
Point(np1) = {lx/2, coord[1],0,m0SE};
Line(10) = {np1, intpListTop[0]};
Line(11) = {np, np1};
Curve Loop(6) = {11,10, 6, 7, -8, -9};
Plane Surface(2) = {6};
Recombine Surface {2};

// Line(9) = {intpListBot[0], intpListTop[0]};
// Curve Loop(10) = {6, 7, -8, 9};
// Plane Surface(2) = {10};

//*** Metal **** // 
//Create box for transfinite
np2 = newp;
coord = Point{intpListTop[nLi + 1]};
Point(np2) = { coord[0], ly/4, 0, m0Li}; 
np3 = newp; 
Point(np3) = {lx/2, ly/4, 0, m0Li};
Line(12) = {np1, np3};
Line(13) = {np3, np2};
Line(14) = {np2, intpListTop[nLi + 1]};
//+
Curve Loop(7) = {12, 13, 14, -6, -10};
//+
Plane Surface(3) = {7};
Recombine Surface {3};//+
Physical Surface("interLayer") = {2};
//+
Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockMetal") = {3};
//+
Physical Curve("blockMetal_top") = {13};
//+
Physical Curve("blockMetal_left") = {14, 7};
//+
Physical Curve("blockMetal_right") = {12, 11};
//+
Physical Curve("blockCeramic_right") = {4};
//+
Physical Curve("blockCeramic_left") = {2};
//+
Physical Curve("blockCeramic_bottom") = {3};
//+
Physical Curve("blockCeramic_top") = {5, 1};
//+
Physical Curve("blockMetal_bottom") = {8, 9};
