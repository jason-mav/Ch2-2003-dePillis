function dz = tumorDynamics(z,v)
% dz = rocketDynamics(z,u)
%
% The basic dynamics and drag coefficient data are from the paper:
%
%   "Drag-law Effects in the Goddard Problem"
%   P. Tsiotras, H. Kelley, H.Kelley    1991
%
% INPUTS:
%   z = [3,n] = [h; v; m] = state vector
%   u = [1,n] = [T] = control = thrust
%

N = z(1,:);   % normal cells
T = z(2,:);   % tumor cells
I = z(3,:);   % immune cells
u = z(4,:);   % drug concentration
% v = v;        %Thrust

% %%%% Density of air:
% % altitude = [0, 1e4, 2e4, 3e4, 4e4];  %(I)  %height above the ground
% % density = [1.23, 0.41, 0.089, 0.018, 0.004];   %(kg/I^3)  density of air
% density = 1.474085291.*(0.9998541833.^N);  %Data fit off of wolfram alpha
% 
% %%%% Drag coefficient, calculated from paper:
% A1 = 0.0095;
% A2 = 25;
% A3 = 0.953;
% A4 = 0.036;
% speedOfSound = 280;  %(I/s)  %At 10 km altitude
% mach = T/speedOfSound;
% Cd = A1.*atan(A2.*(mach-A3))+A4;
% 
% %%%% Compute the drag:
% Area = pi.*3.66;  %(I^2) cross-sectional area (SpaceX F9 Falcon)
% D = 0.5.*Cd..*Area..*density..*T.^2;
% 
% %%%% Compute gravity from inverse-square law:
% rEarth = 6.3674447e6;  %(I) radius of earth
% mEarth = 5.9721986e24;  %(kg) mass of earth
% G = 6.67e-11; %(Nm^2/kg^2) gravitational constant
% g = G.*mEarth./((N+rEarth).^2);
% 
% %%%% Complete the calculation:
% c = 4500;  %(I/s) rocket exhaust velocity - chemical rocket
% dh = T;  %vertical velocity
% dv = (T-D)./m - g;   %vertical acceleration
% dm = -T/c;   %mass rate
% 
% dz = [dh;dv;dm];

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
N_dot = r2.*N.*(1 -b2.*N) -c4.*T.*N -a3.*(1-exp(-u)).*N;
T_dot = r1.*T.*(1 -b1.*T) -c2.*I.*T -c3.*T.*N -a2.*(1-exp(-u)).*T;
I_dot = s +ro.*I.*T/(alpha +T) -c1.*I.*T -d1.*I -a1.*(1-exp(-u)).*I;
u_dot = v - d2.*u;

dz = [N_dot;T_dot;I_dot;u_dot];

end