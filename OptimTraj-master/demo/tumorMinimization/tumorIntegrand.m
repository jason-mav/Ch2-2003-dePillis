function dz = tumorIntegrand(z,v)

% INPUTS:
%   z = [3,n] = [h; v; m] = state vector
%   u = [1,n] = [T] = control = thrust

N = z(1,:);   % normal cells
T = z(2,:);   % tumor cells
I = z(3,:);   % immune cells
u = z(4,:);   % drug concentration
% v = v;        %Thrust

% DePillis 2003 parameters
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

%

% T_dot = r1.*T.*(1 -b1.*T) -c2.*I.*T -c3.*T.*N -a2.*(1-exp(-u)).*T;
T_dot_int = -u.*(T*a2 - T.*r1 + I.*T*c2 + N.*T*c3 + T.^2.*b1.*r1) - T.*a2.*exp(-u); %(int)(T_dot)

dz = T_dot_int;

end