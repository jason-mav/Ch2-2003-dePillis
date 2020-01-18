


%%  Figure 2 !!!
    I = s/d1;
    T = 0.5; 
    N = 0.5;
%     [s,ro] =  meshgrid(0:0.01:35,0:0.01:2);
    ro = linspace(0,2);
    T = 0.5;%linspace(0,1);
    I = 0.8;
    
    s = I.*((c1.*T+d1).*(alpha+T) -ro.*T)/(alpha+T);
%     s = I/(alpha+T).*((c1*T + d1)*(alpha+T) - ro.*T);
    plot(ro,s)
    
    axis([0 2 0 0.35])
    xlabel('?')
    ylabel('s')    
    
%%
%     %% Figure 3
%     clear i I T N xaxis;
%     close all;
%     
%     c1 = 0.78;
%     s = 0.05;
%     
%     step = 0.05; %0.0005;    
%     xaxis = 0:step:2;
%     I(3,length(xaxis)) = 0;
%     T(length(xaxis)) = 0;
%     N(3,length(xaxis)) = 0;
%     
%     for ro=0:step:2
% 
%         syms b;
%         
%         % coexisting equilibrium
%         % (I,T,N) =  (f(b),b,g(b))
%         f_b = s.*(alpha+b)./(c1*b*(alpha+b) + d1*(alpha+b) - ro*b);    
%         g_b = 1 - (c4/r2)*b;
%         
%         % b is a nonnegative solution of: 
%         % b + (c2/r1*b1)*f(b) + (c3/r1*b1)*g(b) - 1/b1 = 0
%         sol = solve(b + (c2/r1*b1)*f_b + (c3/r1*b1)*g_b - 1/b1 == 0, b);    
%         sol_b = vpa(sol);  % has 3 solutions
%         
%         j = 1+int16(ro/step);
%         
%         for i=1:size(sol_b,1)        
%             if sol_b(i)>0 && imag(sol_b(i))==0
%                 
%                 b = real(sol_b(i));        
%                 
%                 I(i,j) = s.*(alpha+b)./(c1*b*(alpha+b) + d1*(alpha+b) - ro*b); % f(b)
%                 T(j)   = b;
%                 N(i,j) = 1 - (c4/r2)*b; % g(b)
%                 
%                 if (I(i,j)>=0 && I(i,j)<=1)
% %                     subplot(3,1,1) 
%                     figure(1)
%                     hold on;
%                     plot(ro,I(i,j),'b.')
%                     xlabel('r')
%                     ylabel('Immune')
%                     hold off;
%                 end
%                 if (T(j)>=0 && T(j)<=1)
% %                     subplot(3,1,2) 
%                     figure(2)
%                     hold on;
%                     plot(ro,T(j),'r.')
%                     xlabel('r')
%                     ylabel('Tumour')
%                     hold off;
%                 end
%                 if (N(i,j)>=0 && N(i,j)<=1)
% %                     subplot(3,1,3) 
%                     figure(3)
%                     hold on;
%                     plot(ro,N(i,j),'b.')
%                     xlabel('r')
%                     ylabel('Normal')
%                     hold off;
%                 end                
%             end
%         end
%        
%     end
    