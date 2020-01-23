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
close all;
save = 1;

%%%% Parameters & Conditions

N_min = 0.75;

v_min = 0;    % Minimum drug dosage
v_max = 1;    % Maximum drug dosage

w1 = 1500;
w2 = 150;
w3 = 1000;
w4 = 40;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Problem Bounds                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initial values
N0 = 1;       % No chemotherapy side effects yet
T0 = 0.25;    % Tumor has already grown

I0 = 0.1001;  % Immune system Low
% I0 = 0.15;    % Immune system High

u0 = 0.01;   % Start the chemo

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
v = soln(end).interp.control(t);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Edit Input!                                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% o.c. might have some negative values
for i=1:length(v)
    if (v(i))<0
        v(i) = 0;
    end
end

max_dose = max(v);
dose_thresh = 0.13*max_dose;
%%2.1*median(v); I0=0.10 
%%795*median(v); %I0 = 0.15
v_bb = v;

% Convert to bang bang
for i=1:length(v)
    if (v(i)) < dose_thresh
        v_bb(i) = 0;
    else
        v_bb(i) = v_max; %max_dose;
    end
end
sprintf('Total drug : %g mg',sum(v_bb))

v_bb_ts = timeseries(v_bb);
I0_ts = timeseries(I0);
total_drug_ts = timeseries(sum(v_bb));

simTime = tf;
sim('model\\model_depillis_bangbang',simTime);



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Print-DirCol!                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

fig1 = figure();
set(gcf,'position',[0 0 700 1000])

subplot(4,1,1);
plot(t,x(1,:), 'LineWidth',1)
axis([0 tf 0 1.5])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Normal cells 10^{11}', 'fontsize',12)
title(sprintf('Minimum normal cells population = %g', min(x(1,:))), 'fontsize',12)

subplot(4,1,2);
plot(t,x(3,:), 'LineWidth',1)
axis([0 tf 0 1.8])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Immune cells 10^{11}', 'fontsize',12)
title(sprintf('Initial immune cells population I_0 = %g', I0), 'fontsize',12)

subplot(4,1,3);
plot(t,x(2,:), 'LineWidth',1)
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Tumor cells 10^{11}', 'fontsize',12)
title(sprintf('Maximum tumor population = %g',  max(x(2,:))), 'fontsize',12)

subplot(4,1,4);
stairs(t,v, 'LineWidth',1)
axis([0 tf 0 1.2])
set(gca, 'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Drug input mg/{L}', 'fontsize',12)
title(sprintf('Total drug : %g mg/L',sum(v)), 'fontsize',12)

% Figure /w hold
fig2 = figure();
hold on;
plot(t,x(1,:),'LineWidth',1)
plot(t,x(2,:), 'LineWidth',1)
plot(t,x(3,:), 'LineWidth',1)
stairs(t,v, 'LineWidth',1)
% plot(t,x(4,:), 'LineWidth',1)
axis([0 tf 0 2])
set(gca,'FontSize',11)
title('Cell Populations and Drug input', 'fontsize',12)
xlabel('Days', 'fontsize',12)
ylabel('Cells (10^{11}), Drug (mg/L)', 'fontsize',12)
legend('N', 'T', 'I', 'v')

% Save the results
if save == 1
    I_0 = int8(I0*100);
    saveas(fig1, sprintf('figures\\I_0=0%d-split', I_0),'fig');
    print(fig1,'-dpng',sprintf('figures\\I_0=0%d-split.png', I_0));

    saveas(fig2, sprintf('figures\\I_0=0%d', I_0),'fig');
    print(fig2,'-dpng',sprintf('figures\\I_0=0%d.png', I_0));
end

sprintf('Total drug given : %g mg',sum(v))
sprintf('Maximum concentration in the body : %g mg/L',max(x(4,:)))


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Print-Bang Bang!                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%%%%% Figure 4x1
fig1_bangbang = figure();
set(gcf,'position',[0 0 700 1000])

subplot(4,1,1);
plot(Cells_out.time,Cells_out.data(:,1), 'LineWidth',1)
axis([0 tf 0 1.5])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Normal cells 10^{11}', 'fontsize',12)
title(sprintf('Minimum normal cells population = %g', min(Cells_out.data(:,1))), 'fontsize',12)

subplot(4,1,2);
plot(Cells_out.time,Cells_out.data(:,3), 'LineWidth',1)
axis([0 tf 0 1.8])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Immune cells 10^{11}', 'fontsize',12)
title(sprintf('Initial immune cells population I_0 = %g', I0), 'fontsize',12)

subplot(4,1,3);
plot(Cells_out.time,Cells_out.data(:,2), 'LineWidth',1)
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Tumor cells 10^{11}', 'fontsize',12)
title(sprintf('Maximum tumor population = %g',  max(Cells_out.data(:,2))), 'fontsize',12)

subplot(4,1,4);
t2 = linspace(1,tf,tf);
stairs(t2, v_bb, 'LineWidth',1)
axis([0 tf 0 1.2])
set(gca, 'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Drug input mg/L', 'fontsize',12)
title(sprintf('Total drug : %g mg/L',sum(v_bb)), 'fontsize',12)

%%%%%%%% Figure /w hold
fig2_bangbang = figure();
hold on;
plot(Cells_out.time,Cells_out.data(:,1), 'LineWidth',1)
plot(Cells_out.time,Cells_out.data(:,2), 'LineWidth',1)
plot(Cells_out.time,Cells_out.data(:,3), 'LineWidth',1)
stairs(t2, v_bb, 'LineWidth',1)
% plot(Cells_out.time,Cells_out.data(:,4), 'LineWidth',1)
axis([0 tf 0 2])
set(gca,'FontSize',11)
title('Bang-Bang approach - Cell Populations and Drug input', 'fontsize',12)
xlabel('Days', 'fontsize',12)
ylabel('Cells (10^{11}), Drug (mg/L)', 'fontsize',12)
legend('N', 'T', 'I', 'v')

% Save the results
if (save == 1)    
    % Print - All populations & Drug input
    I_0 = int8(I0*100);
    saveas(fig1_bangbang, sprintf('figures\\I_0=0%d-split_bangbang', I_0),'fig');
    print(fig1_bangbang,'-dpng',sprintf('figures\\I_0=0%d-split_bangbang.png', I_0));
    
    % Print - All populations & Drug input
    saveas(fig2_bangbang, sprintf('figures\\I_0=0%d_bangbang', I_0),'fig');
    print(fig2_bangbang,'-dpng',sprintf('figures\\I_0=0%d_bangbang.png', I_0));
end

sprintf('Total drug given : %g mg',sum(v))
sprintf('Maximum concentration in the body : %g mg/L',max(x(4,:)))



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Print-Pulsed!                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
fig_pulsed = figure();
hold on;
plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,1), 'LineWidth',1)
plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,2), 'LineWidth',1)
plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,3), 'LineWidth',1)
stairs(pulsed_Cells_out.time, pulsed_Drug_out.data, 'LineWidth',1)
% plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,4), 'LineWidth',1)
axis([0 tf 0 1.2])
set(gca,'FontSize',11)
title('Traditional Pulsed approach - Cell Populations and Drug input', 'fontsize',12)
xlabel('Days', 'fontsize',12)
ylabel('Cells (10^{11}), Drug (mg/L)', 'fontsize',12)
legend('N', 'T', 'I', 'v')

% Save the results
if (save == 1)    
    % Print - All populations & Drug input
    I_0 = int8(I0*100);
    saveas(fig_pulsed, sprintf('figures\\I_0=0%d_pulsed', I_0),'fig');
    print(fig_pulsed,'-dpng',sprintf('figures\\I_0=0%d_pulsed.png', I_0));
    
end

