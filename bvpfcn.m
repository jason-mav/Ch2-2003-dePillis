function p4_dot = bvpfcn(x,p)
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    p4 = p(4);
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    
    p4_dot = -exp(-x4)*(a1*p3*x3 +a2*p2*x2 +a3*p1*x1); % u

%     p_dot = [p2*c3*x1 -p1*(r2 -2*r2*x1^2 -c4*x2 -a3*(1-exp(-x4))) -n; % N    
%             
%             -p3*(ro*alpha*x3/(alpha+x2)^2 -c1*x3) -p2*(r1 -2*r1*b1*x2 -c2*x3 -c3*x1 -a2*(1-exp(-x4))) +p1*c4*x1; % T
%     
%             -p3*(ro*x2/(alpha+x2) -c1*x2 -d1 -a1*(1-exp(-x4))) +p2*c2*x2; % I    
%     
%             -exp(-x4)*(a1*p3*x3 +a2*p2*x2 +a3*p1*x1)]; % u
        
end
        