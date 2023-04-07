lx   = 50;       // width of the model, in unit um
ly   = 50;       // height of the model, in unit um
dw   = 2.0;       // width of the defect located at the top middle, in unit um
dh   = 1.0;       // length of the defect, in unit um
dInt = dw/10;      // Initial thickenss of the interlayer, in unit um
nSE  = 101;      // Number of discretization points of the cosine shape defect for SE
nLi  = 101;      // Number of discretization points of the cosine shape defect for Li
nInt = 101;      // Number of discretization points of the cosine shape defect for interlyer
lx2 = dw*1.1;
lx4 = dw*1.0;

nSE_mesh = 51;
nLi_mesh = 81;

// Create SE geometry


m0SE  = lx/10;    // mesh characteristic length for bottom points of SE
m1SE  = dInt/4;   // mesh characteristic length for interface points of SE
m0Li  = lx/10;    // mesh characteristic length for top points of Li
m1Li  = dInt;   // mesh characteristic length for bottom points of Li
m0Int = dInt/2;   // mesh characteristic length for top points of interlayer
m1Int = dInt/4;   // mesh characteristic length for interface points of interlayer

// *** Bottoom up geometry creation *** //
//** SE geometry */
np_SE_lb = newp;
Point(np_SE_lb) = {0, -ly/2, 0, m0SE};
np_SE_rb = newp;
Point(np_SE_rb) = {lx/2, -ly/2, 0, m0SE};
np_SE_rt = newp;
Point(np_SE_rt) = {lx/2, 0 , 0, m1Int};
For i In {0 : nSE+1}
    x = dw - dw*i/(nSE + 1);
    np= newp;
   SEpList[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Int};
EndFor
Characteristic Length{SEpList[nSE+1]} = m1Int/6;
Line(1) = {np_SE_lb, np_SE_rb};
Line(2) = {np_SE_rb, np_SE_rt};
Line(3) = {np_SE_rt, SEpList[0]};
Spline(4) = SEpList[];
Line(5)  = {SEpList[nSE+1], np_SE_lb};
Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// *** Interlayer **** ///
// ** Interlayer bottom ////
For i In {0 : nSE+1}
    x = dw - dw*i/(nSE + 1);
    np= newp;
   intpListBot[i] = np;
   Point(np)= {x, -dh/2.0*(1+Cos(Pi*x/dw)), 0, m1Int};
EndFor
Characteristic Length {intpListBot[nSE+1]} = m1Int/6;
np_int_rb = newp;
Point(np_int_rb) = {lx2, 0, 0, m1Int};
np_int_rt = newp;
Point(np_int_rt) = {lx2, dInt, 0, m1Int};
np_int_mt = newp;
Point(np_int_mt) = {lx4, dInt, 0, m1Int};
// ** Interlayer top ////
k = 0;
For i In {0 : nSE+1}
    t = i / (nSE + 1);
    x = dw * (1.0 - t);
    y = -dh/2.0 * (1.0 + Cos(Pi * x/dw));
    dx = -dw;
    dy = -dh/2.0 * Sin(Pi * x/dw);
    den = Sqrt(dx*dx + dy*dy);
    newx = x + dInt/den * dy;
    newy = y - dInt/den * dx;
    If (newx >= 0)
        np= newp;
        intpListTop[k] = np;
        Point(np)  = {newx, newy, 0, m1Int};
        k = k + 1;
    EndIf
EndFor
maxTop = k - 1;
Spline(6) = intpListBot[];
Line(7) = {intpListBot[0], np_int_rb};
Line(8) = {np_int_rb, np_int_rt};
Line(9) = {np_int_rt, np_int_mt};
Line(10) = {np_int_mt, intpListTop[0]};
Spline(11) = intpListTop[];
Line(12) = {intpListTop[maxTop], intpListBot[nSE+1]};
Curve Loop(2) = {-6, 7, 8, 9, 10, 11, 12};
Plane Surface(2) = {2};
//+
// **** Li Metal Layer *****
np_li_rt = newp;
Point(np_li_rt) = {lx2, ly/4, 0, m0SE};
np_li_lt = newp;
Point(np_li_lt) = {0, ly/4, 0, m0SE};
np_li_lm = newp;
Point(np_li_lm) = {0, dInt, 0, m1Int};
Line(13) = {np_int_rt, np_li_rt};
Line(14) = {np_li_rt, np_li_lt};
Line(15) = {np_li_lt, np_li_lm};
Line(16) = {np_li_lm, intpListTop[maxTop]};
Curve Loop(3) = {-11, -10, -9, 13, 14, 15, 16};
Plane Surface(3) = {3};
//+
// **** Quartz layer ****
np_q_rt = newp;
Point(np_q_rt) = {lx2, ly/2, 0, m0SE};
np_q_lt = newp;
Point(np_q_lt) = {0, ly/2, 0, m0SE};
Line(17) = {np_li_rt, np_q_rt};
Line(18) = {np_q_rt, np_q_lt};
Line(19) = {np_q_lt, np_li_lt};
Curve Loop(4) = {-14, 17, 18, 19};
Plane Surface(4) = {4};

//+
Physical Surface("interLayer") = {2};
//+
Physical Surface("blockCeramic") = {1};
//+
Physical Surface("blockMetal") = {3};
//+
Physical Surface("blockQuartz") = {4};
// +
Physical Curve("blockMetal_top") = {18};
// +
Physical Curve("blockMetal_bottom") = {6, 7};
// +
Physical Curve("blockCeramic_top") = {3, 4};
// +
Physical Curve("blockCeramic_bottom") = {1};
// +
Physical Curve("blockMetal_left") = {19, 15, 16, 12};
// +
Physical Curve("blockMetal_right") = {8, 13, 17};
// +
Physical Curve("blockCeramic_right") = {2};
// +
Physical Curve("blockCeramic_left") = {5};

Recombine Surface{1, 2, 3, 4};


//+
