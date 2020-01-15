% Figure 2
% coexisting equilibria as a function of ro-s

clear i I T N xaxis;
close all;
% parameters;
% c1 = 1; 

step = 0.01; %0.0005;    
threshold = 1;

for s = 0:step:0.35
% 
%     I(3,length(xaxis)) = 0;
%     T(length(xaxis)) = 0;
%     N(3,length(xaxis)) = 0;

    for ro=0:step*5:2

        syms b f_b g_b;
        
        % coexisting equilibrium
        
        % (I,T,N) =  (f(b),b,g(b))
        
        % f_b = s*(alpha+b)/(c1*b*(alpha+b) + d1*(alpha+b) - ro*b);    
        % g_b = 1 - (c4/r2)*b;

        % b is a nonnegative solution of: 
        % b + (c2/r1*b1)*f(b) + (c3/r1*b1)*g(b) - 1/b1 = 0
        
        f_b = s*(alpha+b)/(c1*b*(alpha+b) + d1*(alpha+b) - ro*b);
        g_b = 1 - (c4/r2)*b;
        
        sol = solve(b + (c2/r1*b1)*f_b + (c3/r1*b1)*g_b - 1/b1 == 0, b);    
        sol_b = vpa(sol);  % has up to 3 solutions

        j=0;  % total number of equilibrium points
        stable_p=0; % number of stable eq. points
        threshold = 1;

        
        for i=1:size(sol_b,1)
            % the solution must be real, positive 
            % and not dividing by zero f_b
            if real(sol_b(i))>0 && imag(sol_b(i))==0 && real(sol_b(i))~=1 
                b = real(sol_b(i));
                f_b = s*(alpha+b)/(c1*b*(alpha+b) + d1*(alpha+b) - ro*b);
                g_b = 1 - (c4/r2)*b;
 
                if b>0 && f_b>0 && g_b>0 % only positive populations make sense
                    
                    j = j+1;

                    b_minus = b - 0.01;
                    f_b_minus = s*(alpha+b_minus)/(c1*b_minus*(alpha+b_minus) + d1*(alpha+b_minus) - ro*b_minus);
                    g_b_minus =  1 - (c4/r2)*b_minus;
                    T = b_minus; I = f_b_minus; N = g_b_minus;
                    I_dot_minus = s +ro*I*T/(alpha +T) -c1*I*T -d1*I;
                    T_dot_minus = r1*T*(1 -b1*T) -c2*I*T -c3*T*N;
                    N_dot_minus = r2*N*(1 -b2*N) -c4*T*N;
                    
                    b_plus  = b + 0.01;
                    f_b_plus = s*(alpha+b_plus)/(c1*b_plus*(alpha+b_plus) + d1*(alpha+b_plus) - ro*b_plus);
                    g_b_plus =  1 - (c4/r2)*b_plus;
                    T = b_plus; I = f_b_plus; N = g_b_plus;
                    I_dot_plus = s +ro*I*T/(alpha +T) -c1*I*T -d1*I;
                    T_dot_plus = r1*T*(1 -b1*T) -c2*I*T -c3*T*N;
                    N_dot_plus = r2*N*(1 -b2*N) -c4*T*N;
                    
                    stable_I = I_dot_minus > 0 && I_dot_plus < 0;
                    stable_T = T_dot_minus > 0 && T_dot_plus < 0;
                    stable_N = N_dot_minus > 0 && N_dot_plus < 0;
                                        
                    if b < 0.4103 %0.4417  %g_b > 0.55  
%                         plot(ro,s,'rx')
%                         hold on;
                          stable_p = stable_p +1;
                    end
                    
% %                     if stable_I
% %                         plot(ro,s,'rx')
% %                         hold on;
% %                     end
% %                     if stable_N
% %                         plot(ro,s,'ko')
% %                         hold on;
% %                     end
% %                     if stable_T
% %                         plot(ro,s,'.','color',[0.8 0.8 0.8])
% %                         hold on;                        
% %                     end
                    
                end
              
            end
        end
        
        unstable_p = j - stable_p;
        
        if j==0 % no equilibrium points ** Region 1 **
            plot(ro,s,'.','color',[0.8 0.8 0.8])
            hold on;
        elseif j==1 % 1 equilibrium point
            hold on;
            if stable_p == 0 % 1 unstable - 0 stable ** Region 2 **
                plot(ro,s,'r.')    
            elseif stable_p == 1 % 1 stable - 0 unstable ** Region 3 **
                plot(ro,s,'ro')
            end
        elseif j==2 % 2 equilibrium points
            hold on;
            if stable_p == 0 % 0 stable - 2 unstable
                plot(ro,s,'b.')
            elseif stable_p == 1 % 1 stable - 1 unstable ** Region 4 **
                plot(ro,s,'bo')
            elseif stable_p == 2 % 2 stable - 0 unstable ** Region 5 **
                plot(ro,s,'ko')
            end
        elseif j==3 % 3 equilibrium points
            hold on;
            if stable_p == 0 % 0 stable - 3 unstable
                plot(ro,s,'g.')
            elseif stable_p == 1 % 1 stable - 2 unstable ** Region 6 **
                plot(ro,s,'go')
            elseif stable_p == 2 % 2 stable - 1 unstable ** Region 7 **
                plot(ro,s,'kx')
            elseif stable_p == 3 % 3 stable - 0 unstable
                plot(ro,s,'g=')
            end
        end

    end
end
axis ([0 2 0 0.35]);
