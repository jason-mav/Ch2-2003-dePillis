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
% Eo = 3.2*10^5;  %cells
% n = 1.101*10^-7;    %1/day/cells

%%
% c2 is larger than the rest
% 0<c3<c2
% Io = s/d1
% To = 10^-5 normalized == 10^6 tumor cells
% No = 1

%% Nullspaces
%
c1 = 1; %0.6;       %0.558 %0.6;
c2 = 0.4875;    %39/80;
c3 = 1.0125;    %81/80;
r1 = 1.5;       %c2+c3
ro = 1; %1;
c4 = 1;

% ok 
alpha = 0.3;
b1 = 1;
b2 = 1; % b2^-1 = 1;
d1 = 0.2;
d2 = 1;
r2 = 1;
s = 0.05;%0.05;

