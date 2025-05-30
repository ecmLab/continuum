lx   = 4;       // width of the model, in unit um
ly   = 4;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.2;       // length of the defect, in unit um
dInt = 0.05;      // Initial thickenss of the interlayer, in unit um
nSE  = 51;      // Number of discretization points of the cosine shape defect for SE
nLi  = 51;      // Number of discretization points of the cosine shape defect for Li
nInt = 51;      // Number of discretization points of the cosine shape defect for interlyer

m0SE  = 1.0;    // mesh characteristic length for bottom points of SE
m1SE  = 0.01;   // mesh characteristic length for interface points of SE
m0Li  = 1.0;    // mesh characteristic length for top points of Li
m1Li  = 0.01;   // mesh characteristic length for bottom points of Li
m0Int = 0.01;   // mesh characteristic length for top points of interlayer
m1Int = 0.01;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
Point(newp) = {0, -ly, 0, m0SE};
Point(newp) = {lx,-ly, 0, m0SE};
bSEp        = newp;     //Start point of the curved boundary of the SE
Point(bSEp) = {lx,  0, 0, m0SE/4};
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

//**** Create Li geometry**/
Point(newp) = {lx, ly, 0, m0Li};
Point(newp) = {0,    ly, 0, m0Li};
bLip        = newp;     //Start point of the curved boundary of the Li
Point(bLip) = {0,  dInt-dh, 0, m1Li};
Line(newl)  = {bLip-2,bLip-1};
Line(newl)  = {bLip-1,bLip};
np          = newp;
nl          = newl;
bLil        = nl;       //Start line of the curved boundary of the Li
For i In {1 : nLi+1}
   x        = dw*i/(nLi + 1);
   Point(np)= {x, dInt-dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Li};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
Point(np)   = {lx, dInt, 0, m1Li};
Line(nl)    = {np-1, np};
eLip        = np;    //End point of the curved boundary of the Li
eLil        = nl;   //End line of the curved boundary of the Li
Line(nl+1)  = {np, bLip-2};
//*** Create loops from Li **/
clLi   = newl;
Curve Loop(clLi) = {bLil-2:eLil+1};

//**** Create Interlayer geometry**/
bIntp       = newp;        //Start from lower-right point of the interlayer
Point(bIntp)= {lx, 0, 0, m1Int*3};
Line(newl)  = {eLip, bIntp};
np          = newp;
nl          = newl;
bIntl       = nl;       //Start line of the lower boundary of the interlayer
For i In {0 : nInt+1}
   x        = dw - dw*i/(nInt + 1);
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Int};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
eIntp       = np-1;    //End at lower-left point of the interlayer
eIntl       = nl-1;   //End line of the lower curved boundary of the interlayer
Line(eIntl+1)= {eIntp, bLip};

//*** Create loops from interlayer **/
clInt   = newl;
Curve Loop(clInt) = {bIntl-1:eIntl+1,bLil:eLil};

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {clSE};
Plane Surface(2) = {clLi};
Plane Surface(3) = {clInt};

//*** Mesh control
// Mesh.RecombineAll = 1;
// Mesh.Algorithm = 6; // Frontal-Delaunay
// Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
// //+
Physical Curve("blockCeramic_bottom") = {bSEl-2};
//+
Physical Curve("blockCeramic_right") = {bSEl-1};
//+
Physical Curve("blockCeramic_top") = {bSEl:eSEl-1};
//+
Physical Curve("blockCeramic_left") = {eSEl};
//+
Physical Curve("blockMetal_top") = {bLil-2};
//+
Physical Curve("blockMetal_left") = {bLil-1,eIntl+1};
//+
Physical Curve("blockMetal_bottom") = {bIntl:eIntl};
//+
Physical Curve("blockMetal_right") = {eLil+1,bIntl-1};
//+
Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockMetal") = {2};
//+
Physical Surface("interLayer") = {3};
//+
Recombine Surface {3, 2,1};
//+
Recombine Surface {3, 2};
//+
Recombine Surface {2};
