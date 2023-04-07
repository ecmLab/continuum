lx   = 0.5;       // width of the model, in unit um
ly   = 0.5;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.2;       // length of the defect, in unit um
dInt = 0.01;      // Initial thickenss of the interlayer, in unit um
nSE  = 51;      // Number of discretization points of the cosine shape defect for SE
nLi  = 51;      // Number of discretization points of the cosine shape defect for Li
nInt = 51;      // Number of discretization points of the cosine shape defect for interlyer

m0SE  = 0.02;    // mesh characteristic length for bottom points of SE
m1SE  = 0.01;   // mesh characteristic length for interface points of SE
m0Li  = 0.05;    // mesh characteristic length for top points of Li
m1Li  = 0.002;   // mesh characteristic length for bottom points of Li
m0Int = 0.005;   // mesh characteristic length for top points of interlayer
m1Int = 0.005;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
Point(newp) = {0, -ly, 0, m0SE};
Point(newp) = {lx,-ly, 0, m0SE};
bSEp        = newp;     //Start point of the curved boundary of the SE
Point(bSEp) = {lx,0, 0, m1SE};
Line(1)     = {bSEp-2, bSEp-1};
Line(2)     = {bSEp-1, bSEp};
For i In {0 : nSE+1}
   x        = dw - dw*i/(nSE + 1);
   np       = newp;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1SE};
EndFor
eSEp        = np;   //End point of the curved boundary of the SE
Line(3)  = {bSEp: eSEp};
Line(4)  = {eSEp, bSEp-2};
//*** Create loops from SE **/
Curve Loop(5) = {1:4};

//**** Create Li geometry**/
Point(newp) = {lx, ly, 0, m0Li};
Point(newp) = {0,    ly, 0, m0Li};
bLip        = newp;     //Start point of the curved boundary of the Li
Point(bLip) = {0,  dInt-dh, 0, m1Li};
Line(11)    = {bLip-2,bLip-1};
Line(12)    = {bLip-1,bLip};
For i In {1 : nLi+1}
   x        = dw*i/(nLi + 1);
   dy       = Pi/2*dh/dw * Sin(Pi*x/dw);
   np       = newp;
   Point(np)= {x-dy*dInt/Sqrt(1+dy*dy), -dh/2.0*(1+Cos(Pi*x/dw))+dInt/Sqrt(1+dy*dy), 0, m1Li};
EndFor
eLip        = newp;     //End point of the curved boundary of the Li
Point(eLip) = {lx, dInt, 0, m1Li};
Line(13)    = {bLip: eLip};
Line(14)    = {eLip, bLip-2};
//*** Create loops from Li **/
Curve Loop(15) = {11:14};

//**** Create Interlayer geometry**/
bIntp       = newp;        //Start from lower-right point of the interlayer
Point(bIntp)= {lx, 0, 0, m1Int};
Line(16)     = {eLip, bIntp};
For i In {0 : nInt+1}
   x        = dw - dw*i/(nInt + 1);
   np       = newp;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Int};
EndFor
eIntp       = np;     //End point of the curved boundary of the interLayer
Line(17) = {bIntp:eIntp};
Line(18) = {eIntp, bLip};
//*** Create loops from interlayer **/
Curve Loop(19) = {16,17,18,13};

//*** Create blocks for the SE, interLayer, and Li-metal **/
Plane Surface(1) = {5};
Plane Surface(3) = {15};
Plane Surface(4) = {19};

//*** Mesh setting
//* Structure mesh on inter layer
Transfinite Surface {4}   = {eIntp,bLip,eLip,bIntp};
Transfinite Curve {13,17} = 151 Using Progression 1;
Transfinite Curve {16,18} = 5 Using Progression 1;
Mesh.RecombineAll = 1;
// Mesh.Algorithm = 8; // Frontal-Delaunay
//Mesh.RecombinationAlgorithm = 3;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockCeramic_bottom") = {1};
//+
Physical Curve("blockCeramic_right") = {2};
//+
Physical Curve("blockCeramic_top") = {3};
//+
Physical Curve("blockCeramic_left") = {4};
//+
Physical Curve("blockMetal_top") = {11};
//+
Physical Curve("blockMetal_left") = {12,18};
//+
Physical Curve("blockMetal_bottom") = {17};
//+
Physical Curve("blockMetal_right") = {14,16};
//+
Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockMetal") = {3};
//+
Physical Surface("interLayer") = {4};
//+

Recombine{4,3,1};//+
Recombine Surface {3};
