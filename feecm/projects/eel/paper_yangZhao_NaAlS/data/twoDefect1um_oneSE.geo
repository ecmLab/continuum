//**** Geometry of the battery model, all in unit mm **/
W = 0.01;       // Width of the battery
La = 0.02;      // thickness of the anode
Le = 0.04;      // thickness of the electrolyte
Lc = 0.02;      // thickness of the cathode
Wd = 0.002;     // width of the cosine shape defect
Ld = 0.001;     // length of the cosine shape defect
Nd  = 101;     // discretization points of the cosine shape defect

e1 = 0.0005;     // coarse mesh control
e2 = 0.00005;     // fine mesh control

//**** Create conner points of each component**/
Point(1) = {0, 0, 0, e1};
Point(2) = {0, W, 0, e1};
Point(3) = {La, W, 0, e2};
Point(4) = {La+Le, W, 0, e1};
Point(5) = {La+Le+Lc, W, 0, e1};
Point(6) = {La+Le+Lc, 0, 0, e1};
Point(7) = {La+Le, 0, 0, e1};
Point(8) = {La, 0, 0, e2};
//**** Create straight lines of each component**/
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {4, 7};
//**** Create points and lines of the first SE/anode defect**/
np          = newp;     //Start point of the defect
nl          = newl;
bSEd        = nl;       //Start line of SE boundary containing the defect
For i In {0 : Nd+1}
   y        = W/4 - Wd/2 + Wd*i/(Nd + 1);
   Point(np)= {La + Ld*Cos(Pi*((y-W/4)/Wd)), y, 0, e2};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor
//**** Create points and lines of the second SE/anode defect**/
For i In {0 : Nd+1}
   y        = 3*W/4 - Wd/2 + Wd*i/(Nd + 1);
   Point(np)= {La + Ld*Cos(Pi*((y-3*W/4)/Wd)), y, 0, e2};
   Line(nl) = {np-1, np};
   np       = newp;
   nl       = newl;
EndFor

eSEd        = nl;   //End line of SE boundary containing the defect
Line(eSEd)  = {np-1, 3};

//*** Create blocks from each component **/
Line Loop(1) = {1, 2, -eSEd:-bSEd, 8};
Line Loop(2) = {bSEd:eSEd, 3, 9, 7};
Line Loop(3) = {-9, 4, 5, 6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

//*** Mesh control
Mesh.CharacteristicLengthMin = 0.00005;
Mesh.CharacteristicLengthMax = 0.003;
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;
Mesh.Algorithm = 8; // Frontal-Quad

//*** Assign names to boundaries and blocks **/
Physical Line("left") = {1};
Physical Line("right") = {5};
Physical Line("top") = {2, 3, 4};
Physical Line("bottom") = {6, 7, 8};
Physical Surface("anode") = {1};
Physical Surface("elyte") = {2};
Physical Surface("cathode") = {3};
