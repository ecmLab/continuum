//**** Create SE geometry**/
Point(1) = {0,   -0.5, 0, 2};
Point(2) = {0.5, -0.5, 0, 2};
Point(3) = {0.5,  0,   0, 2};
Point(4) = {0,    0,   0, 2};
Line(1) = {1, 2};
Line(2) = {2, 3};  
Line(3) = {3, 4};
Line(4) = {4, 1};  
//*** Create loops from SE **/
Curve Loop(1) = {1,2,3,4}; 

//**** Create Li geometry**/
Point(5) = {0,    0.0001, 0, 0.5};
Point(6) = {0.25, 0.0001, 0, 0.5};
Point(7) = {0.25, 0.5,    0, 0.5};
Point(8) = {0,    0.5,    0, 0.5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
//*** Create loops from Li **/
Curve Loop(2) = {5,6,7,8}; 

//*** Create blocks for the SE and Li-metal **/
Plane Surface(1) = {1};
Plane Surface(2) = {2};

//*** Mesh setting
Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Curve {3, 1} = 11 Using Progression 1;
Transfinite Curve {2, 4} = 6 Using Progression 1;
Transfinite Surface {2} = {5, 6, 7, 8};
Transfinite Curve {7, 5} = 61 Using Progression 1;
Transfinite Curve {6, 8} = 21 Using Progression 1;
Recombine Surface {1};
Recombine Surface {2};
//+
//Mesh.RecombineAll = 1;
//Mesh.Algorithm = 8; // Delaunay for quads
//Mesh.RecombinationAlgorithm = 2;

//*** Assign names to boundaries and blocks **/
Physical Curve("blockSE_bottom") = {1};
//+
Physical Curve("blockSE_right")  = {2};
//+
Physical Curve("blockSE_top")    = {3};
//+
Physical Curve("blockSE_left")   = {4};
//+
Physical Curve("blockLi_bottom") = {5};
//+
Physical Curve("blockLi_right")  = {6};
//+
Physical Curve("blockLi_top")    = {7};
//+
Physical Curve("blockLi_left")   = {8};
//+
Physical Surface("blockSE") = {1};
//+
Physical Surface("blockLi") = {2};
