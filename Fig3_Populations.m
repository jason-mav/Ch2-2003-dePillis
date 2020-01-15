    % Figure 3
    clear i I T N xaxis;
    close all;
    parameters;
%     c1=1;
        
    s = 0.33; %0.05;
    
    step = 0.005; %0.0005;
    xaxis = 0:step:2;
    b_step = 0.01;
    I(3,length(xaxis)) = 0;
    T(length(xaxis)) = 0;
    N(3,length(xaxis)) = 0;
    
    for ro=0.01%0:step:2

        syms b;
        
        % coexisting equilibrium
        % (I,T,N) =  (f(b),b,g(b))
        f_b = s.*(alpha+b)./(c1*b*(alpha+b) + d1*(alpha+b) - ro*b);
        g_b = 1 - (c4/r2)*b;
        
        % b is a nonnegative solution of:
        % b + (c2/r1*b1)*f(b) + (c3/r1*b1)*g(b) - 1/b1 = 0
        sol = solve(b + (c2/r1*b1)*f_b + (c3/r1*b1)*g_b - 1/b1 == 0, b);
        sol_b = vpa(sol);  % has 3 solutions
        
        j = 1+int16(ro/step);
        
        for i=1:size(sol_b,1)
            if sol_b(i)>0 %%&& imag(sol_b(i))==0 % solution is real and positive
                                
                b = real(sol_b(i));
                
                N(i,j) = 1 - (c4/r2)*b; % g(b)
                T(j)   = b;
                I(i,j) = s.*(alpha+b)./(c1*b*(alpha+b) + d1*(alpha+b) - ro*b); % f(b)
                
                % will find the stability of the point
                b_minus = b - b_step;
                f_b_minus = s*(alpha+b_minus)/(c1*b_minus*(alpha+b_minus) + d1*(alpha+b_minus) - ro*b_minus);
                g_b_minus =  1 - (c4/r2)*b_minus;
                T_minus = b_minus; I_minus = f_b_minus; N_minus = g_b_minus;
                I_dot_minus = s +ro*I_minus*T_minus/(alpha +T_minus) -c1*I_minus*T_minus -d1*I_minus;
                T_dot_minus = r1*T_minus*(1 -b1*T_minus) -c2*I_minus*T_minus -c3*T_minus*N_minus;
                N_dot_minus = r2*N_minus*(1 -b2*N_minus) -c4*T_minus*N_minus;

                b_plus  = b + b_step;
                f_b_plus = s*(alpha+b_plus)/(c1*b_plus*(alpha+b_plus) + d1*(alpha+b_plus) - ro*b_plus);
                g_b_plus =  1 - (c4/r2)*b_plus;
                T_plus = b_plus; I_plus = f_b_plus; N_plus = g_b_plus;
                I_dot_plus = s +ro*I_plus*T_plus/(alpha +T_plus) -c1*I_plus*T_plus -d1*I_plus;
                T_dot_plus = r1*T_plus*(1 -b1*T_plus) -c2*I_plus*T_plus -c3*T_plus*N_plus;
                N_dot_plus = r2*N_plus*(1 -b2*N_plus) -c4*T_plus*N_plus;

                stable_I = I_dot_minus > 0 && I_dot_plus < 0;
                stable_T = T_dot_minus > 0 && T_dot_plus < 0;
                stable_N = N_dot_minus > 0 && N_dot_plus < 0;

                    figure(1)
                    subplot(3,1,1)                     
                    hold on;
%                     if stable_N
%                         plot(ro,N(i,j),'kx')
%                     else
                        plot(ro,N(i,j),'b.')
%                     end
                    xlabel('r')
                    ylabel('Normal')


%                 if (T(j)>=0 && T(j)<=1)
                    subplot(3,1,2)
                    hold on;
%                     figure(2)
%                     if stable_T
%                        plot(ro,T(j),'kx') 
%                     else
                       plot(ro,T(j),'r.')
%                     end
                    xlabel('r')
                    ylabel('Tumour')
%                 end



%                 if (I(i,j)>=0 && I(i,j)<=1) % cells within limits [0,1]
                    subplot(3,1,3)
                    hold on;
%                     figure(1)
%                     if stable_I
%                         plot(ro,I(i,j),'kx')
%                     else
                        plot(ro,I(i,j),'b.')
%                     end
                    xlabel('r')
                    ylabel('Immune')
                    hold off;
%                 end

                if (T(j)>=0 && T(j)<=1)
                    figure(2)
                    hold on;
                    if stable_T
                       plot(ro,T(j),'kx') 
                    else
                       plot(ro,T(j),'r.')
                    end
                    xlabel('r')
                    ylabel('Tumour')
                    hold off;
                end

            end
        end
       
    end
    