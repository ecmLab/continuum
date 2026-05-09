lx   = 150;       // width of the SE, in unit um
ly   = 150;       // thickness of the SE, in unit um
dw   = 2;       // width of the defect located at the top middle, in unit um
dh   = 2;       // length of the defect, in unit um
nSE  = 101;      // Number of discretization points of the cosine shape defect for SE

m0SE  = 2.0;    // mesh characteristic length for bottom points of SE
m1SE  = 0.05;   // mesh characteristic length for interface points of SE

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
//**** Create the first Ag geometry**/
bAg1p       = np;  //start point of the 1st Ag
bAg1l       = nl;  //the line before the 1st Ag
For i In {0 : nSE+1}
   x        = -dw/4 + dw/2*i/(nSE + 1);
   Point(np)= {x-2*lx/9, dh/2*Cos(Pi*(2*x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eAg1p       = np-1; //end point of the 1st Ag
eAg1l       = nl-1; //end line of the 1st Ag
//**** Create the second Ag geometry**/
bAg2p       = np;  //start point of the 2nd Ag
bAg2l       = nl;  //the line before the 2nd Ag
For i In {0 : nSE+1}
   x        = -dw/4 + dw/2*i/(nSE + 1);
   Point(np)= {x-lx/9, dh/2*Cos(Pi*(2*x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eAg2p       = np-1; //end point of the 2nd Ag
eAg2l       = nl-1; //end line of the 2nd Ag
//**** Create the second defect geometry**/
bD2p        = np;  //start point of the 2nd defect
bD2l        = nl;  //the line before the 2nd defect
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eD2p        = np-1; //end point of the 2nd defect
eD2l        = nl-1; //end line of the 2nd defect
//**** Create the third Ag geometry**/
bAg3p       = np;  //start point of the 3rd Ag
bAg3l       = nl;  //the line before the 3rd Ag
For i In {0 : nSE+1}
   x        = -dw/4 + dw/2*i/(nSE + 1);
   Point(np)= {x+lx/9, dh/2*Cos(Pi*(2*x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eAg3p       = np-1; //end point of the 3rd Ag
eAg3l       = nl-1; //end line of the 3rd Ag
//**** Create the forth Ag geometry**/
bAg4p       = np;  //start point of the 4th Ag
bAg4l       = nl;  //the line before the 4th Ag
For i In {0 : nSE+1}
   x        = -dw/4 + dw/2*i/(nSE + 1);
   Point(np)= {x+2*lx/9, dh/2*Cos(Pi*(2*x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eAg4p       = np-1; //end point of the 4th Ag
eAg4l       = nl-1; //end line of the 4th Ag
//**** Create the third defect geometry**/
bD3p        = np;  //start point of the 3nd defect
bD3l        = nl;  //the line before the 3nd defect
For i In {0 : nSE+1}
   x        = -dw/2 + dw*i/(nSE + 1);
   Point(np)= {x+lx/3, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eD3p        = np-1; //end point of the 3rd defect
eD3l        = nl-1; //end line of the 3rd defect

eSEl        = nl;   //End line of the SE
Line(eSEl)  = {eD3p, 1};
eAg1ll      = newl;   //End line of the 1st Ag
Line(eAg1ll)= {eAg1p, bAg1p};
eAg2ll      = newl;   //End line of the 2nd Ag
Line(eAg2ll)= {eAg2p, bAg2p};
eAg3ll      = newl;   //End line of the 3rd Ag
Line(eAg3ll)= {eAg3p, bAg3p};
eAg4ll      = newl;   //End line of the 4th Ag
Line(eAg4ll)= {eAg4p, bAg4p};

//*** Create loops from SE and Ag **/
clSE        = newl;
Curve Loop(clSE) = {1:eSEl};
clAg1       = newl;
Curve Loop(clAg1) = {bAg1l+1:eAg1l,eAg1ll};
clAg2       = newl;
Curve Loop(clAg2) = {bAg2l+1:eAg2l,eAg2ll};
clAg3       = newl;
Curve Loop(clAg3) = {bAg3l+1:eAg3l,eAg3ll};
clAg4       = newl;
Curve Loop(clAg4) = {bAg4l+1:eAg4l,eAg4ll};

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {clSE};
Plane Surface(2) = {clAg1};
Plane Surface(3) = {clAg2};
Plane Surface(4) = {clAg3};
Plane Surface(5) = {clAg4};

//*** Mesh control
//Mesh.RecombineAll = 1;
//Mesh.Algorithm = 6; // Frontal-Delaunay
//Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockSE_rgt") = {1};
//+
Physical Curve("blockSE_top")   = {2};
//+
Physical Curve("blockSE_lft")  = {3};
//+
Physical Curve("blockSE_btm")   = {bSEl:eD1l+1,eAg1l+1,bD2l:eD2l+1,eAg3l+1,bD3l:eSEl};
//+
Physical Curve("interface")   = {bAg1l+1:eAg1l,bAg2l+1:eAg2l,bAg3l+1:eAg3l,bAg4l+1:eAg4l};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockAg") = {2,3,4,5};

