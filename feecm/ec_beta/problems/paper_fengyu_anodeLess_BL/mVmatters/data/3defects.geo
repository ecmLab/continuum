lx   = 150;       // width of the SE, in unit um
ly   = 150;       // thickness of the SE, in unit um
dw   = 2;       // width of the defect located at the top middle, in unit um
dh   = 2;       // length of the defect, in unit um
nSE  = 101;      // Number of discretization points of the cosine shape defect for SE

m0SE  = 2.0;    // mesh characteristic length for bottom points of SE
m1SE  = 0.2;   // mesh characteristic length for interface points of SE

//**** Create SE geometry**/
Point(1)    = {lx/2, 0,  0, m0SE};
Point(2)    = {lx/2, ly, 0, m0SE};
Point(3)    = {-lx/2,ly, 0, m0SE};
Point(4)    = {-lx/2,0,  0, m0SE};
Line(1)     = {1,2};
Line(2)     = {2,3};
Line(3)     = {3,4};
//**** Create defect geometry**/
//Create the first defect
np          = newp;   
bDfp        = np;
nl          = newl;
bSEl        = nl;       //Start line of SE/CC interface
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x-lx/3, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
//Create the second defect
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
//Create the third defect
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x+lx/3, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eDfp        = np;
eSEl        = nl;   //End line of the SE
Line(eSEl)  = {np-1, 1};

//*** Create loops from SE **/
clSE        = newl;
Curve Loop(clSE) = {1:eSEl};

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {clSE};

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
Physical Curve("blockSE_btm")   = {bSEl:eSEl};
//+
Physical Surface("blockSE") = {1};
