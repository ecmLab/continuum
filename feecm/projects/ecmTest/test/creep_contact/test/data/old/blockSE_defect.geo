lx   = 0.5;       // width of the model, in unit um
ly   = 0.5;       // height of the model, in unit um
dw   = 0.1;    // width of the defect located at the top middle, in unit um
dh   = 0.2;     // length of the defect, in unit um
nSE  = 51;      // Number of discretization points of the cosine shape defect for SE

m0SE  = 0.1;     // mesh characteristic length for bottom points of SE
m1SE  = 0.005;   // mesh characteristic length for interface points of SE

//**** Create SE geometry**/
Point(newp) = {0, -ly, 0, m0SE};
Point(newp) = {lx,-ly, 0, m0SE};
bSEp        = newp;     //Start point of the interface of the SE
Point(bSEp) = {lx,  0, 0, m1SE};
Line(newl)  = {bSEp-2, bSEp-1};
Line(newl)  = {bSEp-1, bSEp};  
np          = newp;
nl          = newl;
bSEl        = nl;       //Start line of the interface of the SE
For i In {0 : nSE+1}
   x        = dw - dw*i/(nSE + 1);
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
Line(nl)  = {np-1, bSEp-2};
//*** Create loops from SE **/
clSE        = newl;
Curve Loop(clSE) = {bSEl-2:nl}; 

//*** Create blocks for the SE **/
Plane Surface(1) = {clSE};

//*** Mesh setting
Mesh.RecombineAll = 1;
Mesh.Algorithm = 8; // Delaunay for quads
Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockSE_bottom") = {bSEl-2};
//+
Physical Curve("blockSE_right") = {bSEl-1};
//+
Physical Curve("blockSE_left") = {nl};
//+
Physical Curve("blockSE_top") = {bSEl:nl-1};
//+
Physical Surface("blockSE") = {1};
