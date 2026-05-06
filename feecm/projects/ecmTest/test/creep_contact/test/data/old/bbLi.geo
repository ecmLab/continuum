//**** Create Li geometry**/
Point(5) = {0,    0.0001, 0, 0.5};
Point(6) = {0.5, 0.0001, 0, 0.5};
Point(7) = {0.5, 0.5,    0, 0.5};
Point(8) = {0,    0.5,    0, 0.5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
//*** Create loops from Li **/
Curve Loop(2) = {5,6,7,8}; 

//*** Create blocks for Li-metal **/
Plane Surface(2) = {2};

//*** Mesh setting
Transfinite Surface {2} = {5, 6, 7, 8};
Transfinite Curve {7, 5} = 61 Using Progression 1;
Transfinite Curve {6, 8} = 21 Using Progression 1;
Recombine Surface {2};

//*** Assign names to boundaries and blocks **/
//Physical Curve("blockLi_bottom") = {5};
//+
//Physical Curve("blockLi_right")  = {6};
//+
//Physical Curve("blockLi_top")    = {7};
//+
//Physical Curve("blockLi_left")   = {8};
//+
//Physical Surface("blockLi") = {2};
