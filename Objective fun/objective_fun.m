% function tumor = objective_fun(x)

%   Copyright 2004 The MathWorks, Inc.  

    
    
%   tumor - objective function
fun = @(x) -x(4)*(x(2)*a2 - x(2)*r1 + x(3)*x(2)*c2 + x(1)*x(2)*c3 + x(2)^2*b1*r1) -x(2)*a2*exp(-x(4)); %int(r1*x(2)*(1-b1*x(2)) -c2*x(1)*x(2) -c3*x(2)*x(3) -a2*(1-exp(-x(4)))*x(2));
    

% %   constraints, ceq(x)=0
% 
% %   N_dot = f_N , N =(int)f_N , -N +(int)f_N = 0
% g1 = @(x) -x(1) -x(4)*(x(1)*a3 - x(1)*r2 + x(1)*x(2)*c4 + x(1)^2*b2*r2) -x(1)*a3*exp(-x(4)); 
% 
% %   (x(1) - 0.75);
% %     g2 = @(x) x(2);
%     
% %   -I +(int)f_I = 0
% g3 = @(x) -x(3) -x(3)*a1*exp(-x(4)) -(x(4)*(x(3)*x(2)*a1 - alpha*s - x(2)*s + x(3)*x(2)*d1 + x(3)*a1*alpha + x(3)*alpha*d1 - x(3)*x(2)*ro + x(3)*x(2)^2*c1 + x(3)*x(2)*alpha*c1))/(x(2) + alpha);
%     
% %   -u +(int)f_u = 0
% g4 = @(x) -x(4) +(v*(v - 2*d2*x(4)))/2;
%     
%     fx = @(x) f(x(:,1),x(:,2),x(:,3)) + 
%          1E20 * (x(:,1))

x0 = [1;0.2;0.15;0];
A = [];
b = [];
Aeq = []; 
beq = [];    

% lower & upper bounds
lb = [0.75;0;0;0];
ub = [inf;inf;inf;1]; % v has the bound, not x4

nlcon = @nonlcon;

x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nlcon)