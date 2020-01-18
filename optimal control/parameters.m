%% DePillis 2003
a1 = 0.2;
a2 = 0.3;
a3 = 0.1;

b1 = 1;     %0.3030; %1;
b2 = 1;     % b2^-1 = 1;

alpha = 0.3;

c1 = 1;     %0.08
c2 = 0.5;   %1667;  %0.5 %39/80;
c3 = 1;     %81/80;
c4 = 1;     %0.3;   % 1

d1 = 0.2;   %0.066;  %0.2
d2 = 1;

r1 = 1.5;   %c2+c3
r2 = 1;

s = 0.33;   %0.33 %0.05;

ro = 0.01;  %1;



% Nondiamensionalization
% b1 = alpha*b1;
% 
% c1 = alpha*c1/r2;
% c2 = c2*s/r2;
% c3 = c3/b2/r2;
% c4 = alpha*c4/r2;
% 
% d1 = d1/r2;
% 
% r1 = r1/r2;
% 
% ro = ro/r2;



%% Kuznetzov
% 
% r1 = 1.636;     %1/day - ?
% b1 = 2*10^-3;   %1/cells - ? 
% s = 0.1181;     %cells/day - ?
% ro = 1.131;     %1/day - ?
% alpha = 20.19;  %cells - ?
% c1 = 0.00311;   %1/day/cells - ?
% d1 = 0.3743;    %1/day - ?
% 
% n = 1.101*10^-7;    %1/day/cells

%% DePillis 2001
% c2 is larger than the rest
% 0<c3<c2
% Io = s/d1
% To = 10^-5 (normalized *10^11) == 10^6 tumor cells
% No = 1

