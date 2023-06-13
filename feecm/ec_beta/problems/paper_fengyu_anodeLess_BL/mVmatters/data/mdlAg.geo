lx   = 40;       // width of the SE, in unit um
ly   = 40;       // thickness of the SE, in unit um
dw   = 2;       // width of the defect located at the top middle, in unit um
dh   = 2;       // length of the defect, in unit um
nSE  = 101;      // Number of discretization points of the cosine shape defect for SE

m0SE  = 0.5;    // mesh characteristic length for bottom points of SE
m1SE  = 0.05;   // mesh characteristic length for interface points of SE

//**** Create SE geometry**/
Point(1)    = {lx/2, 0,  0, m0SE};
Point(2)    = {lx/2, ly, 0, m0SE};
Point(3)    = {-lx/2,ly, 0, m0SE};
Point(4) = {-lx/2,0,  0, m0SE};
Line(1)     = {1,2};
Line(2)     = {2,3};
Line(3)     = {3,4};
np          = newp;     //Start point of SE bottom boundary is most left point
nl          = newl;
bSEl        = nl;       //Start line of the curved boundary of the SE
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eSEl        = nl;   //End line of the SE
Line(eSEl)  = {np-1, 1};
//*** Create loops for SE **/
clSE        = newl;
Curve Loop(clSE) = {1:nl};

//*** Create loops for Ag particle **/
eAgl        = newl;   //Start line of the Ag
Line(eAgl)  = {np-1, 5};
clAg        = newl;
Curve Loop(clAg) = {5:nl-1, eAgl};

//*** Create blocks for the SE and Ag particle **/
Plane Surface(1) = {clSE};
Plane Surface(2) = {clAg};

//*** Mesh control
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockSE_rgt") = {1};
//+
Physical Curve("blockSE_top")   = {2};
//+
Physical Curve("blockSE_lft")  = {3};
//+
Physical Curve("blockSE_btm_lft")   = {4};
//+
Physical Curve("interface")     = {5:nl-1};
//+
Physical Curve("blockSE_btm_rgt")   = {eSEl};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockAg") = {2};
