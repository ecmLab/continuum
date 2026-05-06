//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1.0, 0, 0, 1.0};
//+
Point(3) = {1.0, 0.58, 0, 1.0};
//+
Point(4) = {0.0, -1.0, 0, 1.0};
//+
Point(5) = {1.0, -1.0, 0, 1.0};
//+
Point(6) = {0.0, 0.1, 0, 1.0};
//+
Point(7) = {1.0, 0.68, 0, 1.0};
//+
Point(8) = {0, 1.0, 0, 1.0};
//+
Point(9) = {1.0, 1.0, 0, 1.0};
//+
Point(10) = {0, 0.0, 0, 1.0};
//+
Point(11) = {1.0, 0.58, 0, 1.0};

//+
Line(1) = {4, 5};
//+
Line(2) = {5, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 4};
//+
Line(5) = {10, 11};
//+
Line(6) = {11, 7};
//+
Line(7) = {7, 6};
//+
Line(8) = {6, 10};

//+
Line(9) = {6, 8};
//+
Line(10) = {8, 9};
//+
Line(11) = {9, 7};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {8, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {9, 10, 11, 7};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {4, 2, 1, 3} = 2 Using Progression 1;
//+
Transfinite Curve {8, 6, 11, 9} = 2 Using Progression 1;
//+
Transfinite Curve {10, 7, 5} = 4 Using Progression 1;
