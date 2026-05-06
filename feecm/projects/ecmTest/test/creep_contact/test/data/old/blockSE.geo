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

//*** Create blocks for the SE**/
Plane Surface(1) = {1};

//*** Mesh setting
Transfinite Surface {1} = {1, 2, 3, 4};
Transfinite Curve {3, 1} = 11 Using Progression 1;
Transfinite Curve {2, 4} = 6 Using Progression 1;
Recombine Surface {1};

