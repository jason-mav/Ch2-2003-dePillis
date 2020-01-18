syms s alpha b1 c1 c2 c3 c4 d1 r1 r2 ro I T N f_T;

% sol = solve ( I == s*(alpha+T)/(c1*T*(alpha+T) + d1*(alpha+T) - ro*T), c1)

alpha = 0.3;
b1 = 1;
b2 = 1; % b2^-1 = 1;
d1 = 0.2;
d2 = 1;
r2 = 1;

s = 0.05;%0.05;
ro = 1; % ??
% r1 = 2; % ??
% c4 = 1;


% solve(10*(13*c1)/10 == 37/50, c1)
% a = diff(s.*(alpha+T)./(c1.*T.*(alpha+T) + d1.*(alpha+T) - ro.*T),T)
% T = 1
% d1 = 0.2;
% alpha = 0.3; 
% solve( (13*c1 - 37/5) )

%% Figure 1
% 
% %% N1
% % (T,I) = 
% % (1, 2.5)
% % (0.85, -4)
% % (0.1, -1)
% % (0.5, 3)
% % (0, 0.8)
% 
% %% Approximate c1 (= 0.6)
% T=1; I=2.5;
% sol = solve ( I == s*(alpha+T)/(c1*T*(alpha+T) + d1*(alpha+T) - ro*T), c1)
% T=0.85; I=-4;
% sol = solve ( I == s*(alpha+T)/(c1*T*(alpha+T) + d1*(alpha+T) - ro*T), c1)
% T=0.1; I=-1;
% sol = solve ( I == s*(alpha+T)/(c1*T*(alpha+T) + d1*(alpha+T) - ro*T), c1)
% T=0.5; I=3;
% sol = solve ( I == s*(alpha+T)/(c1*T*(alpha+T) + d1*(alpha+T) - ro*T), c1)
% T=0; I=0.8;
% sol = solve ( I == s*(alpha+T)/(c1*T*(alpha+T) + d1*(alpha+T) - ro*T), c1)
% 
% %% N2-N3
% % T = 1/b1 - (c2/r1*b1).*I - (c3/r1*b1).*N;
% % N = 1 - (c4/r2).*T;
% 
% %% N2 points
% I=0; T=-0.35; N=2;
% sol = solve ( T == 1/b1 - (c2/r1*b1)*I - (c3/r1*b1)*N, c3) % c3= (27*r1)/40
% 
% I=2; T=-1; N=2;
% sol = solve ( T == 1/b1 - (c2/r1*b1)*I - (c3/r1*b1)*N, c3) % c3= r1 -c2
% 
% I=2; T=0; N=0;
% sol = solve ( T == 1/b1 - (c2/r1*b1)*I - (c3/r1*b1)*N, c2) % c2= r1/2
% 
% I=0; T=0.5; N=0.8;
% sol = solve ( T == 1/b1 - (c2/r1*b1)*I - (c3/r1*b1)*N, c3) 
% 
% I=1; T=0.001; N=1;
% sol = solve ( T == 1/b1 - (c2/r1*b1)*I - (c3/r1*b1)*N, c3) % c2=(55*r1)/148 -(145*c3)/148
% 
% 
% eqns = [ c3 == (27*r1)/40, c3 == r1 -c2, c3 == 5*r1/8];
% vars = [ c2 c3 r1];
% [c2, c3, r1] = solve(eqns, vars)
% 
% %% N3 points
% I=0; T=2; N=-1;
% sol = solve ( N == 1 - (c4/r2)*T, c4)
% 
% I=2; T=2; N=-1;
% sol = solve ( N == 1 - (c4/r2)*T, c4)
% 
% I=2; T=0; N=1;
% sol = solve ( N == 1 - (c4/r2)*T, c4)
% 
% I=1; T=0.001; N=1;
% sol = solve ( N == 1 - (c4/r2)*T, c4)

%%








