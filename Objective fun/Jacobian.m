syms x y z c1 c2 c3 c4 e ro d r

% 2003 p.16
% x_dot = x*(1 -x) -c4*x*y -f1*x                     = -x*(f1 + x + c4*y - 1)
% y_dot = r*y*(1 -b*y) -c2*y*z -c3*x*y -f2*y         = -y*(f2 - r + c3*x + c2*z + b*r*y)
% z_dot = 1 + ro*y*z/(1+y) - c1*y*z -d*z -f3*z       = (ro*y*z)/(y + 1) - f3*z - c1*y*z - d*z + 1

% f1 = a1*(1 -exp(-u))
% f2 = a2*(1 -exp(-u))
% f3 = a3*(1 -exp(-u))

% No drug - Unperturbed
F = [ x*(1 -x)-c4*x*y, r*y*(1 -b*y)-c2*y*z-c3*x*y, 1+(ro*y*z)/(1+y)-c1*y*z-d*z ];
v = [ x y z ];
J_u = jacobian(F,v);

% simplified 
J_u = [ 1-c4*y-2*x,   -c4*x,                  0;
      -c3*y,        r-2*b*r*y-c2*z-c3*x,    -c2*y;
      0,            (ro*z)/((1+y)^2)-c1*z,   (ro*y)/(1+y)-c1*y-d];


  
% With drug - Perturbed
F = [ x*(1 -x)-c4*x*y-e*x, r*y*(1 -b*y)-c2*y*z-c3*x*y-e*y, 1+(ro*y*z)/(1+y)-c1*y*z-d*z-e*z ];
v = [ x y z ];
J_p = jacobian(F,v);
  
% simplified 
J_p = [ 1-c4*y-2*x-e,   -c4*x,                      0;
        -c3*y,          r-2*b*r*y-c2*z-c3*x-e,      -c2*y;
        0,              (ro*z)/((1+y)^2)-c1*z,      (ro*y)/(1+y)-c1*y-d-e];

  