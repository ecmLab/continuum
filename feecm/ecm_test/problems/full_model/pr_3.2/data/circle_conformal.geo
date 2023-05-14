lx   = 0.5;       // width of the model, in unit um
ly   = 0.5;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.3;       // length of the defect, in unit um
dInt = 0.01;      // Initial thickenss of the interlayer, in unit um

m0SE  = 0.05;    // mesh characteristic length for bottom points of SE
m1SE  = 0.005;   // mesh characteristic length for interface points of SE
m0Li  = 0.01;    // mesh characteristic length for top points of Li
m1Li  = 0.002;   // mesh characteristic length for bottom points of Li
m0Int = 0.005;   // mesh characteristic length for top points of interlayer
m1Int = 0.001;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
Point(1) = {0,    -ly, 0, m0SE};
Point(2) = {lx,   -ly, 0, m0SE};
Point(3) = {lx,     0, 0, m1SE};
Point(4) = {dw,     0, 0, m1SE};
Point(5) = {dw, dw-dh, 0, m1SE};
Point(6) = {0,  dw-dh, 0, m1SE};
Point(7) = {0,    -dh, 0, m1SE};
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 5};
Circle(5)  = {5, 6, 7};
Line(6)  = {7, 1};
//*** Create loops from SE **/
Curve Loop(7) = {1:6};

//**** Create Li geometry**/
Point(8) = {lx/2,    ly, 0, m0Li};
Point(9) = {0,       ly, 0, m0Li};
Point(10) = {0,  dInt-dh, 0, m1Li};
Point(11) = {0,  dw-dh, 0, m1Li};
Point(12) = {dw-dInt,  dw-dh, 0, m1Li};
Point(13) = {dw-dInt,    dInt, 0, m1Li};
Point(14) = {lx/2,  dInt, 0, m1Li};
Line(8)  = {8, 9};
Line(9)  = {9, 10};
Circle(10)={10,11,12};
Line(11)  = {12, 13};
Line(12) = {13, 14};
Line(13) = {14, 8};
//*** Create loops from Li **/
Curve Loop(14) = {8:13};

//**** Create Interlayer geometry**/
Point(15)= {lx/2,0,  0, m1Int}; //Start from lower-right point of the interlayer
Point(16)= {dw,  0,  0, m1Int}; 
Point(17)= {dw, dw-dh,  0, m1Int}; 
Point(18)= {0,  dw-dh,  0, m1Int};
Point(19)= {0,    -dh,  0, m1Int};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 17};
Circle(18)={17,18,19};
Line(19) = {19, 10};

//*** Create loops from interlayer **/
Curve Loop(20) = {15:19,10:12};

//*** Create blocks for the SE interlayer, and Li-metal **/
Plane Surface(1) = {7};
Plane Surface(2) = {14};
Plane Surface(3) = {20};

//*** Mesh control
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
//+
Physical Curve("blockCeramic_bottom") = {1};
//+
Physical Curve("blockCeramic_right") = {2};
//+
Physical Curve("blockCeramic_top") = {3:5};
//+
Physical Curve("blockCeramic_left") = {6};
//+
Physical Curve("blockMetal_top") = {8};
//+
Physical Curve("blockMetal_left") = {9,19};
//+
Physical Curve("blockMetal_bottom") = {16:18};
//+
Physical Curve("blockMetal_right") = {13,15};
//+
Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockMetal") = {2};
//+
Physical Surface("interLayer") = {3};
//+
Recombine Surface {3, 2, 1};
//+
Recombine Surface {3, 2, 1};
//+
Recombine Surface {2};
