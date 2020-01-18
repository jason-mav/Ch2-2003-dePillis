    
parameters;
syms b;

% coexisting equilibriua
% (I,T,N) =  (f(b),b,g(b))
f_b = s.*(alpha+b)./(c1*b*(alpha+b) + d1*(alpha+b) - ro*b);    
g_b = 1 - (c4/r2)*b;

% b is a nonnegative solution of: 
% b + (c2/r1*b1)*f(b) + (c3/r1*b1)*g(b) - 1/b1 = 0
sol = solve(b + (c2/r1*b1)*f_b + (c3/r1*b1)*g_b - 1/b1 == 0, b);    
sol_b = vpa(sol);  % has 3 solutions

for i=1:size(sol_b,1)        
    if sol_b(i)>0 && imag(sol_b(i))==0

        b = real(sol_b(i));
        f_b = s.*(alpha+b)./(c1*b*(alpha+b) + d1*(alpha+b) - ro*b); % f(b)
        g_b = 1 - (c4/r2)*b; % g(b)
        fprintf('Coexisting Equilibrium at (N,T,I)=(%g,%g,%g)',g_b,b,f_b)
        if b > 0.4103
            fprintf(' - stable.\n')
        else
            fprintf(' - unstable.\n')
        end
    end
end

% dead equilibriua
syms a;
% (I,T,N) =  (f(a),a,g(a))
f_a = s.*(alpha+a)./(c1*a*(alpha+a) + d1*(alpha+a) - ro*a);    
g_a = 1 - (c4/r2)*a;

% a is a nonnegative solution of: 
% a + (c2/r1*b1)*f(a) + (c3/r1*b1)*g(a) - 1/b1 = 0
sol = solve(a + (c2/r1*b1)*f_a - 1/b1 == 0, a);    
sol_a = vpa(sol);  % has 3 solutions

fprintf('Dead Equilibrium at (N,T,I)=(%g,%g,%g) - Type 1\n',0,0,s/d1)

for i=1:size(sol_a,1)        
    if sol_a(i)>0 && imag(sol_a(i))==0

        a = real(sol_a(i));
        f_a = s.*(alpha+a)./(c1*a*(alpha+a) + d1*(alpha+a) - ro*a); % f(a)
        g_a = 1 - (c4/r2)*a; % g(a)
        fprintf('Dead Equilibrium at (N,T,I)=(%g,%g,%g) - Type 2\n',g_a,a,f_a)
        
    end
end


