lx   = 0.5;       // width of the model, in unit um
ly   = 0.5;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.2;       // length of the defect, in unit um
dInt = 0.01;      // Initial thickenss of the interlayer, in unit um
nSE  = 51;      // Number of discretization points of the cosine shape defect for SE
nLi  = 51;      // Number of discretization points of the cosine shape defect for Li
nInt = 51;      // Number of discretization points of the cosine shape defect for interlyer

m0SE  = 0.2;    // mesh characteristic length for bottom points of SE
m1SE  = 0.05;   // mesh characteristic length for interface points of SE
m0Li  = 0.02;    // mesh characteristic length for top points of Li
m1Li  = 0.001;   // mesh characteristic length for bottom points of Li
m0Int = 0.02;   // mesh characteristic length for top points of interlayer
m1Int = 0.01;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
// Point(newp) = {0, -ly, 0, m0SE};
// Point(newp) = {lx,-ly, 0, m0SE};
// bSEp        = newp;     //Start point of the curved boundary of the SE
// Point(bSEp) = {lx,  0, 0, m1SE};
// Line(newl)  = {bSEp-2, bSEp-1};
// Line(newl)  = {bSEp-1, bSEp};
// np          = newp;
// nl          = newl;
// bSEl        = nl;       //Start line of the curved boundary of the SE
// For i In {0 : nSE+1}
//    x        = dw - dw*i/(nSE + 1);
//    Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
//    Line(nl) = {np-1, np};
//    np       = newp;
//    nl       = newl;
// EndFor
// eSEl        = nl;   //End line of the SE
// Line(eSEl)  = {np-1, bSEp-2};
// //*** Create loops from SE **/
// clSE        = newl;
// Curve Loop(clSE) = {bSEl-2:nl};

//**** Create Li geometry**/
Point(newp) = {lx/2, ly, 0, m0Li};
Point(newp) = {0,    ly, 0, m0Li};
bLip        = newp;     //Start point of the curved boundary of the Li
Point(bLip) = {0,  dInt-dh, 0, m1Li};
Line(1)  = {bLip-2,bLip-1};
Line(2)  = {bLip-1,bLip};
For i In {1 : nLi+1}
   x        = dw*i/(nLi + 1);
   dy       = Pi/2*dh/dw * Sin(Pi*x/dw);
   np       = newp;
   Point(np)= {x-dy*dInt/Sqrt(1+dy*dy), -dh/2.0*(1+Cos(Pi*x/dw))+dInt/Sqrt(1+dy*dy), 0, m1Li};
EndFor
eLip        = newp;     //End point of the curved boundary of the Li
Point(eLip) = {lx/2, dInt, 0, m1Li};
Line(3)  = {bLip: eLip};
Line(4)  = {eLip, bLip-2};
//*** Create loops from Li **/
Curve Loop(5) = {1:4};

//**** Create Interlayer geometry**/
bIntp       = newp;        //Start from lower-right point of the interlayer
Point(bIntp)= {lx/2, 0, 0, m1Int};
Line(6)     = {eLip, bIntp};
For i In {0 : nInt+1}
   x        = dw - dw*i/(nInt + 1);
   np       = newp;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Int};
EndFor
eIntp       = np;     //End point of the curved boundary of the interLayer
Line(7) = {bIntp:eIntp};
Line(8) = {eIntp, bLip};

//*** Create loops from interlayer **/
Curve Loop(9) = {6,7,8,3};

//*** Create blocks for the SE and Li-metal **/
// Plane Surface(1) = {clSE};
Plane Surface(2) = {5};
Plane Surface(3) = {9};

//*** Mesh setting
// Transfinite Surface {3} = {eIntp,bLip,eLip,bIntp};
// Transfinite Curve {3} = 201 Using Progression 1;
// Transfinite Curve {7} = 201 Using Progression 1;
// Transfinite Curve {6,8} = 5 Using Progression 1;
Mesh.RecombineAll = 1;
//Mesh.Algorithm = 8; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;
// Recombine Surface {3,2};

//*** Assign names to boundaries and blocks **/
//+
//Physical Curve("blockCeramic_bottom") = {bSEl-2};
//+
//Physical Curve("blockCeramic_right") = {bSEl-1};
//+
//Physical Curve("blockCeramic_top") = {bSEl:eSEl-1};
//+
//Physical Curve("blockCeramic_left") = {eSEl};
//+
Physical Curve("blockMetal_top") = {1};
//+
Physical Curve("blockMetal_left") = {2,8};
//+
Physical Curve("blockMetal_bottom") = {7};
//+
Physical Curve("blockMetal_right") = {4,6};
//+
// Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockMetal") = {2};
//+
Physical Surface("interLayer") = {3};
