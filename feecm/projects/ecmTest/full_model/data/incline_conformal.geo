lx   = 0.5;       // width of the model, in unit um
ly   = 0.5;       // height of the model, in unit um
dw   = 0.1;       // width of the defect located at the top middle, in unit um
dh   = 0.1;       // length of the defect, in unit um
dInt = 0.01;      // Initial thickenss of the interlayer, in unit um

m0SE  = 0.05;    // mesh characteristic length for bottom points of SE
m1SE  = 0.01;   // mesh characteristic length for interface points of SE
m0Li  = 0.02;    // mesh characteristic length for top points of Li
m1Li  = 0.005;   // mesh characteristic length for bottom points of Li
m0Int = 0.01;   // mesh characteristic length for top points of interlayer
m1Int = 0.005;   // mesh characteristic length for interface points of interlayer

//**** Create SE geometry**/
Point(1) = {0, -ly, 0, m0SE};
Point(2) = {lx,-ly, 0, m0SE};
Point(3) = {lx,  0, 0, m1SE};
Point(4) = {dw,  0, 0, m1SE};
Point(5) = {0, -dh, 0, m1SE};
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 5};
Line(5)  = {5, 1};
//*** Create loops from SE **/
Curve Loop(6) = {1:5};

//**** Create Li geometry**/
Point(6) = {lx/2,    ly, 0, m0Li};
Point(7) = {0,       ly, 0, m0Li};
Point(8) = {0,  dInt-dh, 0, m1Li};
Point(9) = {dw,    dInt, 0, m1Li};
Point(10)= {lx/2,  dInt, 0, m1Li};
Line(7)  = {6, 7};
Line(8)  = {7, 8};
Line(9)  = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 6};
//*** Create loops from Li **/
Curve Loop(12) = {7:11};

//**** Create Interlayer geometry**/
Point(11)= {lx/2,0,  0, m1Int}; //Start from lower-right point of the interlayer
Point(12)= {dw,  0,  0, m1Int}; 
Point(13)= {0,  -dh, 0, m1Int}; 
Line(13) = {10, 11};
Line(14) = {11, 12};
Line(15) = {12, 13};
Line(16) = {13, 8};

//*** Create loops from interlayer **/
Curve Loop(17) = {13:16,9,10};

//*** Create blocks for the SE interlayer, and Li-metal **/
Plane Surface(1) = {6};
Plane Surface(2) = {12};
Plane Surface(3) = {17};

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
Physical Curve("blockCeramic_top") = {3,4};
//+
Physical Curve("blockCeramic_left") = {5};
//+
Physical Curve("blockMetal_top") = {7};
//+
Physical Curve("blockMetal_left") = {8,16};
//+
Physical Curve("blockMetal_bottom") = {14,15};
//+
Physical Curve("blockMetal_right") = {11,13};
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
