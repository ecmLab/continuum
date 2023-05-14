// ** define size and mesh quality
lx	= 0.005;	//width of the LPS, in unit mm
ly	= 0.005;	//height of the LPS, in unit mm
rc	= 0.0025;	//radius of the NMC, in unit mm
m1      = 0.0001;  //mesh size for out boundaries
m2      = 0.00002; //mesh size for interface

//** create NMC geometry **/
Point(1) 	= {0,	0,	0,	m1};
Point(2) 	= {rc,	0,	0,	m2};
Point(3) 	= {0,	rc,	0,	m2};
Line(1)	 	= {1,2};
Circle(2)	= {2,1,3};
Line(3)		= {3,1};
//** create loop from NMC **/
Curve Loop(4) = {1:3};

//** create LPS geometry **/
Point(4) 	= {rc,	0,	0,	m2};
Point(5) 	= {lx,	0,	0,	m1};
Point(6) 	= {lx,	ly,	0,	m1};
Point(7) 	= {0,	ly,	0,	m1};
Point(8) 	= {0,	rc,	0,	m2};
Point(9)	= {0,	0,	0,	m1};
Line(5)		= {4,5};
Line(6)		= {5,6};
Line(7)		= {6,7};
Line(8)		= {7,8};
Circle(9)	= {8,9,4};

//** create loop from LPS **/
Curve Loop(10) = {5:8,9};


//*** Create blocks for the NMC and LPS **/
Plane Surface(1) = {4};
Plane Surface(2) = {10};

//*** Mesh control
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;

//** Assign names to boundaries and blocks **/
//+
Physical Curve("block_NMC_right") = {2};
//+
Physical Curve("block_LPS_left") = {9};
//+
Physical Curve("block_bottom") = {1,5};
//+
Physical Curve("block_right") = {6};
//+
Physical Curve("block_top") = {7};
//+
Physical Curve("block_left") = {8,3};
//+
Physical Surface("block_NMC") = {1};
//+
Physical Surface("block_LPS") = {2};
//+

