%% Function to find common tangent of two polynomial functions
% p1: parameters of function 1; p2: parameters of function 2
% iniVec: initial vectors
function [P1, P2] = comTangent(p1, p2, iniVec)
    syms x eqm1 eqm2 x1 x2;
    f1(x)    = p1(1) + p1(2)*x + p1(3)*x^2 + p1(4)*x^3 + p1(5)*x^4 + p1(6)*x^5;
    f2(x)    = p2(1) + p2(2)*x + p2(3)*x^2 + p2(4)*x^3 + p2(5)*x^4 + p2(6)*x^5;
    f1u(x)   = diff(f1(x),x);
    f2u(x)   = diff(f2(x),x);
    eqm1     = (f2(x2)-f1(x1))/(x2-x1) - f1u(x1) == 0;
    eqm2     = f2u(x2) - f1u(x1) == 0;
    [tmp1,tmp2] = vpasolve([eqm1 eqm2], [x1 x2], iniVec);
    P1(:,1)  = double(tmp1(real(tmp1)>0&imag(tmp1)==0));
    P1(:,2)  = double(subs(f1(x),x,P1(:,1)));
    P2(:,1)  = double(tmp2(real(tmp2)>0&imag(tmp2)==0));
    P2(:,2)  = double(subs(f2(x),x,P2(:,1)));
    