//** Representative Volume Element (RVE) for Cathode Expansion with Coating Material **//


lx	= 0.0050;	//width of the LPS, in unit mm
ly	= 0.0050;	//height of the LPS, in unit mm
th 	= 0.000080;  //thickness of the coating material, in unit mm
rc	= 0.0025;	//radius of the NMC, in unit mm


//** create NMC geometry **/
Point(1) 	= {0,	0,	0,	0.1};
Point(2) 	= {rc,	0,	0,	0.05};
Point(3) 	= {0,	rc,	0,	0.05};
Line(1)	 	= {1,2};
Circle(2)	= {2,1,3};
Line(3)		= {3,1};

//** create loop from NMC **/
Curve Loop(4) = {1:3};


//** create LPS geometry **/
Point(4) 	= {rc+th,	0,	0,	0.05};
Point(5) 	= {lx,	0,	0,	0.1};
Point(6) 	= {lx,	ly,	0,	0.1};
Point(7) 	= {0,	ly,	0,	0.1};
Point(8) 	= {0,	rc+th,	0,	0.05};
Point(9)	= {0,	0,	0,	0.1};
Line(5)		= {4,5};
Line(6)		= {5,6};
Line(7)		= {6,7};
Line(8)		= {7,8};
Circle(9)	= {8,9,4};

//** create loop from LPS **/
Curve Loop(10) = {5:9};


//** create coating material geometry **/
Point(10) 	= {0,	rc,	0,	0.05};
Point(11) 	= {rc,	0,	0,	0.05};
Point(12)	= {0,	0,	0,	0.1};
Line(11)	= {11,4};
Line(12)	= {8,10};
Circle(13)	= {10,12,11};

//** create loop from coating material **/
Curve Loop(14) = {11,-9,12,13};


//*** Create blocks for the NMC and LPS **/
Plane Surface(1) = {4};
Plane Surface(2) = {10};
Plane Surface(3) = {14};


//*** Mesh control
Mesh.RecombineAll = 1;
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2;

//** Assign names to boundaries and blocks **/
//+
Physical Curve("block_master") = {2};
//+
Physical Curve("block_slave") = {13};
//+
Physical Curve("block_bottom") = {1,11,5};
//+
Physical Curve("block_right") = {6};
//+
Physical Curve("block_top") = {7};
//+
Physical Curve("block_left") = {8,12,3};
//+
Physical Surface("block_NMC") = {1};
//+
Physical Surface("block_LPS") = {2};
//+
Physical Surface("block_coating") = {3};
//+


