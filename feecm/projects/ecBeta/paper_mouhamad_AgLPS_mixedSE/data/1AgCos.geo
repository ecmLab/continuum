lx   = 100;       // width of the SE, in unit um
ly   = 100;       // thickness of the SE, in unit um
dw   = 2;        // width of the pre-existing crack located at the bottom middle, in unit um
dh   = 10;       // length of the pre-existing crack, in unit um
nCk  = 100;      // Number of discretization points of the cosine shape defect for SE

m0SE  = 2.0;    // mesh characteristic length for corner points of SE
m1SE  = 0.01;   // mesh characteristic length for defect points of SE

//**** Create SE geometry**/
Point(1)    = {dw/2, 0, 0, m1SE};
Point(2)    = {lx/2, 0,  0, m0SE};
Point(3)    = {lx/2, ly, 0, m0SE};
Point(4)    = {0,    ly, 0, m0SE};
Point(5)    = {0,    dh, 0, m1SE};
Line(1)     = {1,2};
Line(2)     = {2,3};
Line(3)     = {3,4};
Line(4)     = {4,5};
//**** Create the pre-existing crack geometry**/
bCp         = 5;   //start point of the crack
bCl         = 5;   //start line of the crack
np          = newp;
nl          = newl;
For i In {1 : nCk-1}
   x        = dw/2*i/nCk;
   Point(np)= {x, dh*Cos(Pi*(x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eCp         = np-1;      //end point of the crack
Line(nl)    = {eCp,1};   // the last segment
eCl         = nl;        //end line of the crack

//*** Create loops for SE and Ag **/
clSE        = newl;
Curve Loop(clSE) = {1:nl};

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {clSE};

//*** Mesh control
//Mesh.RecombineAll = 1;
//Mesh.Algorithm = 6; // Frontal-Delaunay
//Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
//Physical Curve("blockSE_btm") = {1};
//+
Physical Curve("blockSE_rgt") = {2};
//+
Physical Curve("blockSE_top") = {3};
//+
Physical Curve("blockSE_lft") = {4};
//+
Physical Curve("liMetal")   = {1,eCl-nCk/10:99,101:eCl};
//+
Physical Curve("ag")   = {100};
//+
Physical Surface("blockSE") = {1};

