//**** Create SE geometry**/
dx = 9.8;
dm = 0.5;
Point(1) = {-dx/2,  -dx/2, 1.5*(dx-0.8), dm};
Point(2) = {dx/2,   -dx/2, 1.5*(dx-0.8), dm};
Point(3) = {dx/2,    dx/2, 1.5*(dx-0.8), dm};
Point(4) = {-dx/2,   dx/2, 1.5*(dx-0.8), dm};
Line(1) = {1, 2};
Line(2) = {2, 3};  
Line(3) = {3, 4};
Line(4) = {4, 1};  
//*** Create loops from lines **/
Curve Loop(1) = {1,2,3,4}; 

//*** Create blocks for the insert surface **/
Plane Surface(1) = {1};

//*** Assign names to boundaries and blocks **/
Physical Curve("bottom") = {1};
//+
Physical Curve("right")  = {2};
//+
Physical Curve("top")    = {3};
//
Physical Curve("left")   = {4};
//+
Physical Surface("Ins")  = {1};
