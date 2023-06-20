lx   = 100;       // width of the SE, in unit um
ly   = 100;       // thickness of the SE, in unit um
dw   = 2;       // width of the defect located at the top middle, in unit um
dh   = 2;       // length of the defect, in unit um
nSE  = 101;      // Number of discretization points of the cosine shape defect for SE

m0SE  = 1.0;    // mesh characteristic length for bottom points of SE
m1SE  = 0.1;   // mesh characteristic length for interface points of SE

//**** Create SE geometry**/
Point(1)    = {lx/2, 0,  0, m0SE};
Point(2)    = {lx/2, ly, 0, m0SE};
Point(3)    = {-lx/2,ly, 0, m0SE};
Point(4)    = {-lx/2,0,  0, m0SE};
Line(1)     = {1,2};
Line(2)     = {2,3};
Line(3)     = {3,4};
//**** Create the first defect geometry**/
np          = newp;
bD1p        = np;   //start point of 1st defect
nl          = newl;
eD1l        = nl;   //the line before the 1st defect
bSEl        = nl;   //start line of the CC/SE interface
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x-lx/3, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eD1p        = np-1; //end point of the 1st defect
eD1l        = nl-1; //end line of the 1st defect
//**** Create Ag geometry**/
//Create the first Ag particle
bAgp        = np;  //start point of the 1st Ag
bAgl        = nl;  //the line before the 1st Ag
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eAgp        = np-1; //end point of the 1st Ag
eAgl        = nl-1; //end line of the 1st Ag
//**** Create the second defect geometry**/
bD2p        = np;  //start point of the 2nd defect
bD2l        = nl;  //the line before the 2nd defect
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x+lx/3, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eD2p        = np-1; //end point of the 2nd defect
eD2l        = nl-1; //end line of the 2nd defect

eSEl        = nl;   //End line of the SE
Line(eSEl)  = {eD2p, 1};
eAgll       = newl;   //End line of the Ag
Line(eAgll) = {eAgp, bAgp};

//*** Create loops from SE and Ag **/
clSE        = newl;
Curve Loop(clSE) = {1:eSEl};
clAg        = newl;
Curve Loop(clAg) = {bAgl+1:eAgl,eAgll};

//*** Create blocks for the SE and Li-metal **/
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
Physical Curve("blockSE_btm_lft")   = {bSEl:bAgl};
//+
Physical Curve("blockSE_btm_rgt")   = {bD2l:eSEl};
//+
Physical Curve("interface")   = {bAgl+1:eAgl};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockAg") = {2};

