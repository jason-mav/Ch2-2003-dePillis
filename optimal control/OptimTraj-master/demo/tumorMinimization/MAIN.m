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

%%%% Parameters & Conditions

N_min = 0.75;

v_min = 0;    % Minimum drug dosage
v_max = 1;    % Maximum drug dosage

w1 = 1500;
w2 = 150;
w3 = 1000;
w4 = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initial values
N0 = 1;       % No chemotherapy side effects yet
T0 = 0.25;    % Tumor has already grown

I0 = 0.1001;  % Immune system Low
% I0 = 0.15;    % Immune system High

u0 = v_max;   % Start the chemo

% Final desired values
Nf = 1;     % Healthy after treatment
Tf = 0;     % Trying to eradicate the tumor
If = 1.65;  % Immune population of a healthy organism
uf = 0;     % End of treatment

tf = 150;   % Duration of chemo (days)

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
v_Low = v_min; % Minimum dosage
v_Upp = v_max; % Maximum dosage

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

P.bounds.control.low = v_Low;
P.bounds.control.upp = v_Upp;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Initial Guess                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% guess the progression of the tumor
P.guess.time = [0, tf];  %(s)
P.guess.state = [ [N0;T0;I0;u0],  [Nf;Tf;If;uf] ];
P.guess.control = [v_Upp, v_Low];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Objective and Dynamic functions                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Dynamics function:
P.func.dynamics = @(t,x,u)( tumorDynamics(x,u) );

% Objective function:
% P.func.bndObj = @(t0,x0,tF,xF)( xF(2) );  % Minimize tumor // Maximize final height -xF(1)/10000

P.func.pathObj = @(t,x,u)( w1*x(2,end) +w2*tumorIntegrand(x) +w3*max(x(2,:)) +w4*u );
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
        P.options.method = 'trapezoid';
        P.options.defaultAccuracy = 'high';
        P.options.nlpOpt.MaxFunEvals = 2e5;
        P.options.nlpOpt.MaxIter = 1e5;
        P.options.trapezoid.nGrid = 50;
        
    case 'rungeKutta'
        P.options.method = 'rungeKutta';
        P.options.defaultAccuracy = 'high';
        P.options.nlpOpt.MaxFunEvals = 5e4;
        P.options.nlpOpt.MaxIter = 2e2;

    case 'chebyshev'            
        P.options.method = 'chebyshev';
        P.options.defaultAccuracy = 'low';
        P.options.chebyshev.nColPts = 15;
    
    case 'hermiteSimpson'        
        P.options.method = 'hermiteSimpson';
        P.options.defaultAccuracy = 'high';        
        P.options.hermiteSimpson.nSegment = 50;        
        P.options.nlpOpt.MaxFunEvals = 5e4;
        P.options.nlpOpt.MaxIter = 2e2;
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

% o.c. might have some negative values
for i=1:length(u)
    if (u(i))<0
        u(i) = 0;
    end
end
drug_thresh = sum(u);

fig1 = figure();
set(gcf,'position',[0 0 700 1000])

subplot(4,1,1);
plot(t,x(1,:), 'LineWidth',1)
axis([0 tf 0 1.5])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Normal cells', 'fontsize',12)
title(sprintf('Minimum normal cells population = %g', min(x(1,:))), 'fontsize',12)

subplot(4,1,2);
plot(t,x(3,:), 'LineWidth',1)
axis([0 tf 0 1.8])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Immune cells', 'fontsize',12)
title(sprintf('Io = %g', I0), 'fontsize',12)

subplot(4,1,3);
plot(t,x(2,:), 'LineWidth',1)
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Tumor cells', 'fontsize',12)
title(sprintf('Maximum tumor population = %g',  max(x(2,:))), 'fontsize',12)

subplot(4,1,4);
stairs(t,u, 'LineWidth',1)
axis([0 tf 0 1.2])
set(gca, 'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Drug input', 'fontsize',12)
title(sprintf('Total drug : %g ?mg/mL',sum(u)), 'fontsize',12)

I_0 = int8(I0*100);
saveas(fig1, sprintf('figures\\I_0=0%d-split', I_0),'fig');
print(fig1,'-dpng',sprintf('figures\\I_0=0%d-split.png', I_0));

fig2 = figure();
hold on;
plot(t,x(1,:),'LineWidth',1)
plot(t,x(2,:), 'LineWidth',1)
plot(t,x(3,:), 'LineWidth',1)
stairs(t,u, 'LineWidth',1)
% plot(t,x(4,:), 'LineWidth',1)
axis([0 tf 0 2])
set(gca,'FontSize',11)
title('Cell Populations and Drug input', 'fontsize',12)
xlabel('Days', 'fontsize',12)
ylabel('Cells', 'fontsize',12)
legend('N', 'T', 'I', 'v')

sprintf('Total drug given : %g \t??mg/mL??',sum(u))
sprintf('Maximum concentration in the body : %g \t??mg/mL??',max(x(4,:)))

saveas(fig2, sprintf('figures\\I_0=0%d', I_0),'fig');
print(fig2,'-dpng',sprintf('figures\\I_0=0%d.png', I_0));

simTime = tf;
I_0 = double(I0);

% normal_cells = x(1,:);
% tumor_cells = x(2,:);
% immune_cells = x(3,:);
% drug_conc = x(4,:);
% drug_input = u;




