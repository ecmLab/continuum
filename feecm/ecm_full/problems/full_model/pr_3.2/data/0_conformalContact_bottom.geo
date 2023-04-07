lx   = 0.5;       // width of the model, in unit um
ly   = 0.5;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.2;       // length of the defect, in unit um
dInt = 0.01;      // Initial thickenss of the interlayer, in unit um
nSE  = 51;      // Number of discretization points of the cosine shape defect for SE
nLi  = 51;      // Number of discretization points of the cosine shape defect for Li
nInt = 51;      // Number of discretization points of the cosine shape defect for interlyer

m0SE  = 0.05;    // mesh characteristic length for bottom points of SE
m1SE  = 0.005;   // mesh characteristic length for interface points of SE
m0Li  = 0.02;    // mesh characteristic length for top points of Li
m1Li  = 0.002;   // mesh characteristic length for bottom points of Li
m0Int = 0.002;   // mesh characteristic length for top points of interlayer
m1Int = 0.001;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
Point(newp) = {0, -ly, 0, m0SE};
Point(newp) = {lx,-ly, 0, m0SE};
bSEp        = newp;     //Start point of the curved boundary of the SE
Point(bSEp) = {lx,  0, 0, m1SE};
Line(newl)  = {bSEp-2, bSEp-1};
Line(newl)  = {bSEp-1, bSEp};
np          = newp;
nl          = newl;
bSEl        = nl;       //Start line of the curved boundary of the SE
For i In {0 : nSE+1}
   x        = dw - dw*i/(nSE + 1);
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eSEl        = nl;   //End line of the SE
Line(eSEl)  = {np-1, bSEp-2};
//*** Create loops from SE **/
clSE        = newl;
Curve Loop(clSE) = {bSEl-2:nl};

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {clSE};


//*** Mesh control
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockCeramic_bottom") = {bSEl-2};
//+
Physical Curve("blockCeramic_right") = {bSEl-1};
//+
Physical Curve("blockCeramic_left") = {eSEl};
//+
Physical Curve("blockCeramic_top") = {bSEl:eSEl-1};

Recombine Surface {1};
