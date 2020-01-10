% MAIN.m -- Goddard Rocket
%
% This script runs a trajectory optimization to find the optimal thrust
% trajectory for the rocket to reach the maximum altitude. Physical
% parameters are roughly based on the SpaceX Falcon 9 rocket.
%
% Dynamics include variable mass, inverse-square gravity, speed-dependent
% drag coefficient, height dependent air density.
%
% NOTES:
%   This problem sort of converges, but not very well. I think that there
%   is a singular arc in it that is not being handled correctly. It is
%   still interesting to see as an example of ways in which problems might
%   misbehave.
%

clc; clear;
addpath ../../

%%%% Assumptions of parameters
parameters;
N_min = 0.75;
v_min = 0;    % Minimum drug dosage
v_max = 1;    % Maximum drug dosage
w1 = 1500;
w2 = 150;
w3 = 1000;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initial values
N0 = 1;     %Rocket starts on the ground
T0 = 0.2;   %Rocket starts stationary
I0 = 0.15;  %Rocket starts full of fuel
u0 = v_max;     %Rocket starts full of fuel

% Final desired values
Tf = 0;     %Trying to eradicate the tumor
tf = 150;

% mF = mEmpty; %Assume that we use all of the fuel

% Cells population

% Normal
N_Low = N_min;
N_Upp = inf; 

% Tumor
T_Low = 0;    %Just look at the trajectory as it goes up
T_Upp = inf;  % Go as fast as you can

% Immune
I_Low = 0;
I_Upp = inf;

% drug concentration
u_Low = 0;
u_Upp = inf; % practically == 1 

% drug input
vLow = v_min; % Minimum dosage
vUpp = v_max; % Maximum dosage

P.bounds.initialTime.low = 0;
P.bounds.initialTime.upp = 0;

P.bounds.finalTime.low = tf;
P.bounds.finalTime.upp = tf;

P.bounds.state.low = [N_Low;T_Low;I_Low;u_Low];
P.bounds.state.upp = [N_Upp;T_Upp;I_Upp;u_Upp];

P.bounds.initialState.low = [N0;T0;I0;u0];
P.bounds.initialState.upp = [N0;T0;I0;u0];

P.bounds.finalState.low = [N_Low;Tf;I_Low;u_Low];
P.bounds.finalState.upp = [N_Upp;Tf;I_Upp;u_Upp];

P.bounds.control.low = vLow;
P.bounds.control.upp = vUpp;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Initial Guess                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% hGuess = 2e4;   %(m) guess at the maximum height reached
P.guess.time = [0, tf];  %(s)
P.guess.state = [ [N0;T0;I0;u0],  [1;Tf;1.6;0] ];
P.guess.control = [vUpp, vLow];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Objective and Dynamic functions                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Dynamics function:
P.func.dynamics = @(t,x,u)( tumorDynamics(x,u) );

% Objective function:
% P.func.bndObj = @(t0,x0,tF,xF)( xF(2) );  % Minimize tumor // Maximize final height -xF(1)/10000

P.func.pathObj = @(t,x,u)( w1*x(2,end) +w2*tumorIntegrand(x,u) +w3*max(x(2,:)) +u );
% P.func.pathObj = @(t,x,u)( ones(size(t)) ); % minimize time

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Options and Method selection                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% method = 'trapezoid';
% method = 'rungeKutta';
% method = 'chebyshev';
method = 'hermiteSimpson';

switch method
    
    case 'trapezoid'
        
%         P.options(1).method = 'trapezoid';
%         P.options(1).defaultAccuracy = 'low';
        
        P.options.method = 'trapezoid';
        P.options.defaultAccuracy = 'high';
        P.options.nlpOpt.MaxFunEvals = 2e5;
        P.options.nlpOpt.MaxIter = 1e5;
        P.options.trapezoid.nGrid = 50;
        
    case 'rungeKutta'
        P.options(1).method = 'rungeKutta';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'rungeKutta';
        P.options(2).defaultAccuracy = 'medium';
        
    case 'chebyshev'
        
        P.options(1).method = 'chebyshev';
        P.options(1).defaultAccuracy = 'low';
        
        P.options(2).method = 'chebyshev';
        P.options(2).defaultAccuracy = 'low';
        P.options(2).chebyshev.nColPts = 15;
    
    case 'hermiteSimpson'        
        P.options.method = 'hermiteSimpson';
        P.options.defaultAccuracy = 'high';
        
        P.options.nlpOpt.MaxFunEvals = 5e4;
        P.options.nlpOpt.MaxIter = 2e4;
        P.options.hermiteSimpson.nSegment = 50;
end


%%%% NOTES:
%
% 1) Orthogonal collocation (chebyshev) is not a good method for this problem, beause there is a
% discontinuity in solution of the thrust curve. It still sort of works,
% but will find a sub-optimal answer, or produce ringing.
%
% 2) Why does the 'trapezoid' low resolution version finish so quickly and the medium
% quality one take forever? Hint: Look at the feasibility printout: it is
% cyclical. If you were to plot the solution as a function of iteration,
% you would find that occasionally the discontinuity moves, which causes a
% consistency error in the NLP. Eventually it gets to the "right" answer,
% although it is pretty boring. I suspect that you could get more
% interesting behavior with different constants.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Solve!                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln = optimTraj(P);

t = linspace(soln(end).grid.time(1),soln(end).grid.time(end),tf);
x = soln(end).interp.state(t);
u = soln(end).interp.control(t);

figure();
subplot(2,2,1);
plot(t,x(1,:))
axis([0 tf 0 1.5])
xlabel('time (s)')
ylabel('Normal cells')
title(sprintf('minimum %g',min(x(1,:))))
% title('Maximal Height Trajectory')
subplot(2,2,2);
plot(t,x(3,:))
axis([0 tf 0 1.8])
xlabel('time (s)')
ylabel('Immune cells')
% title('Goddard Rocket')
subplot(2,2,3);
plot(t,x(2,:))
xlabel('time (s)')
ylabel('Tumor cells')
subplot(2,2,4);
plot(t,u)
axis([0 tf 0 1.2])
xlabel('time (s)')
ylabel('Drug input')
title(sprintf('Total drug : %g ?mg/mL',sum(x(4,:))))

% figure()
% hold on;
% plot(t,x(1,:))
% plot(t,x(2,:))
% plot(t,x(3,:))
% plot(t,x(4,:))
% plot(t,u)
% legend('Normal cells', 'Tumor cells', 'Immune cells', 'Drug concentration', 'Drug input')

sprintf('Total drug given : %g \t??mg/mL??',sum(x(4,:)))
sprintf('Maximum concentration : %g \t??mg/mL??',max(x(4,:)))

% plot(t,x(4,:))
% xlabel('time (s)')
% ylabel('Drug concentration')
