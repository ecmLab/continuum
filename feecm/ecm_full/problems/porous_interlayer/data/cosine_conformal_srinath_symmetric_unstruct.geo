lx   = 50;       // width of the model, in unit um
ly   = 50;       // height of the model, in unit um
dw   = 2.0;       // width of the defect located at the top middle, in unit um
dh   = 6.0;       // length of the defect, in unit um
dInt = dw/10;      // Initial thickenss of the interlayer, in unit um
nSE  = 101;      // Number of discretization points of the cosine shape defect for SE
nLi  = 101;      // Number of discretization points of the cosine shape defect for Li
nInt = 101;      // Number of discretization points of the cosine shape defect for interlyer

nSE_mesh = 51; 
nLi_mesh = 81;

// Create SE geometry


m0SE  = lx/10;    // mesh characteristic length for bottom points of SE
m1SE  = dInt/2;   // mesh characteristic length for interface points of SE
m0Li  = lx/10;    // mesh characteristic length for top points of Li
m1Li  = dInt;   // mesh characteristic length for bottom points of Li
m0Int = dInt/2;   // mesh characteristic length for top points of interlayer
m1Int = dInt/3;   // mesh characteristic length for interface points of interlayer

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
k = 0;
For i In {0 : nLi+1}
    t = i / (nLi + 1);
    x = dw * (1.0 - t);
    y = -dh/2.0 * (1.0 + Cos(Pi * x/dw));
    dx = -dw;
    dy = -dh/2.0 * Sin(Pi * x/dw);
    den = Sqrt(dx*dx + dy*dy);
    newx = x + dInt/den * dy;
    newy = y - dInt/den * dx;
    If (newx >= 0)
        np= newp;
        intpListTop[k] = np;
        Point(np)  = {newx, newy, 0, m1Int};
        k = k + 1;
    EndIf
EndFor

// np = newp;
// intpListTop[k] = np;
// Point(np) = {0, -dh + dInt, 0, m0Int};
Spline(6) = intpListTop[];

For i In {0 : nLi+1}
    x = dw - dw*i/(nLi + 1);
    np= newp;
   intpListBot[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Int};
   // Point(np)= {x, 0, 0, m0Int};
EndFor
maxTop = k-1;

Line(7) = {intpListTop[maxTop], intpListBot[nLi+1]};
Spline(8) = intpListBot[];
np = newp;
Point(np) = {lx/2, 0, 0, m0Int};
npx = newp; 
Point(npx) = {lx/10, 0, 0, m0Int};
Line(901) = {np, npx};
Line(9) = {npx, intpListBot[0]};

np1 = newp;
coord = Point{intpListTop[0]};
Point(np1) = {lx/2, coord[1],0,m0Int};
np2 = newp; 
Point(np2) = {lx/10, coord[1],0, m0Int};
Line(101) = {np1, np2};
Line(10) = {np2, intpListTop[0]};
Line(11) = {np, np1};
Curve Loop(6) = {11,101, 10, 6, 7, -8, -901, -9};
Plane Surface(2) = {6};
Recombine Surface {2};

//*** Metal **** // 
//Create box for transfinite
np2 = newp;
coord = Point{intpListTop[maxTop]};
Point(np2) = { coord[0], ly/4, 0, m0Li}; 
np3 = newp; 
Point(np3) = {lx/2, ly/4, 0, m0Li};
Line(12) = {np1, np3};
Line(13) = {np3, np2};
np3 = newp;
coord = Point{intpListTop[maxTop]};
Point(np3) = {0, 0, 0, m1Int*4};
Line(14) = {np2, np3};
Line(15) = {np3, intpListTop[maxTop]};
//+
Curve Loop(7) = {12, 13, 14,15, -6, -101, -10};
//+
Plane Surface(3) = {7};
//+
Recombine Surface {3};
//+
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
Physical Curve("blockMetal_bottom") = {8, 901, 9};
//+
Characteristic Length {59} = 0.025;
//+
Characteristic Length {209} = 0.01;
