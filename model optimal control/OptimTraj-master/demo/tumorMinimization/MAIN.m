% MAIN.m
% 
% Create a model for the cancerous area in the patient's body and simulate
% chemotherapy treatment based on DIRCOL, bang bang and traditionally pulsed regimens.
% 

clc; clear;
addpath ../../
close all;

save = 0; % Save the printed figures


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

I0 = 0.10;  % Immune system Low
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
method = 'hermiteSimpson';

switch method
    
    case 'trapezoid'        
        P.options.method = 'trapezoid';
        P.options.defaultAccuracy = 'high';
        P.options.nlpOpt.MaxFunEvals = 2e5;
        P.options.nlpOpt.MaxIter = 1e5;
        P.options.trapezoid.nGrid = 50;
        
    case 'hermiteSimpson'        
        P.options.method = 'hermiteSimpson';
        P.options.defaultAccuracy = 'high';        
        P.options.hermiteSimpson.nSegment = 149;        
%         P.options.nlpOpt.MaxFunEvals = 5e4;
        P.options.nlpOpt.MaxIter = 50; %2e2;
end


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
if I0 == 0.15
    dose_thresh = 0.12*max_dose; 
else % 0.10
    dose_thresh = 0.1455*max_dose; 
end
close all;

% Convert to bang bang
v_bb = v;
for i=1:length(v)
    if (v(i)) < dose_thresh
        v_bb(i) = 0;
    else
        v_bb(i) = v_max; %max_dose;
    end
end

% Convert to timeseries
v_bb_ts = timeseries(v_bb);
I0_ts = timeseries(I0);
total_drug_ts = timeseries(sum(v_bb));

simTime = tf;
sim('model\\model_depillis_bangbang',simTime);


fprintf('[DirCol] Total drug given : %g mg/m^2 \n',sum(v)) 
fprintf('[DirCol] Maximum concentration in the body : %g mg/L \n',max(x(4,:)))

fprintf('[Bang Bang] Dose threshold: %g mg/m^2\n', dose_thresh)
fprintf('[Bang Bang] Total drug given : %g mg/m^2 \n',sum(v_bb))
fprintf('[Bang Bang] Maximum concentration in the body : %g mg/L \n',max(Cells_out.data(:,4)))

% I0=0.10 : 
% [DirCol] Total drug given : 15.8379 mg/m^2 
% [DirCol] Maximum concentration in the body : 0.722654 mg/L 
% [Bang Bang] Dose threshold: 0.145499 mg/m^2
% [Bang Bang] Total drug given : 16 mg/m^2 
% [Bang Bang] Maximum concentration in the body : 0.98603 mg/L 

% I0=0.15 : 
% [DirCol] Total drug given : 14.9764 mg/m^2 
% [DirCol] Maximum concentration in the body : 0.660461 mg/L 
% [Bang Bang] Dose threshold: 0.119961 mg/m^2
% [Bang Bang] Total drug given : 15 mg/m^2 
% [Bang Bang] Maximum concentration in the body : 0.997847 mg/L 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Print-DirCol!                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
I_0 = int8(I0*100);

fig1 = figure();
set(gcf,'position',[0 0 700 1000])

subplot(4,1,1);
plot(t,x(1,:), 'LineWidth',1)
axis([0 tf 0 1.5])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Normal cells (10^{11})', 'fontsize',12)
title(sprintf('Minimum normal cells population = %g', min(x(1,:))), 'fontsize',12)

subplot(4,1,2);
plot(t,x(3,:), 'LineWidth',1)
axis([0 tf 0 1.8])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Immune cells (10^{11})', 'fontsize',12)
title(sprintf('Initial immune cells population I_0 = 0.%d', I_0), 'fontsize',12)

subplot(4,1,3);
plot(t,x(2,:), 'LineWidth',1)
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Tumor cells (10^{11})', 'fontsize',12)
title(sprintf('Maximum tumor cells population = %g',  max(x(2,:))), 'fontsize',12)

subplot(4,1,4);
plot(t,x(4,:), 'LineWidth',1)
axis([0 tf 0 1.2])
set(gca, 'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Drug concentration (mg/L)', 'fontsize',12)
title(sprintf('Maximum drug concentration: %g',max(x(4,:))), 'fontsize',12)

% hold on;
% stairs(t,v, 'LineWidth',1)
% axis([0 tf 0 1.2])
% set(gca, 'FontSize',11)
% xlabel('Days', 'fontsize',12)
% ylabel('Drug input (mg/m^2)', 'fontsize',12)
% title(sprintf('Total drug : %g mg/L',sum(v)), 'fontsize',12)

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
ylabel('Cells (10^{11}), Drug (mg/m^2)', 'fontsize',12)
legend('N', 'T', 'I', 'v')

% Save the results
if save == 1    
    saveas(fig1, sprintf('figures\\I_0=0%d-split', I_0),'fig');
    print(fig1,'-dpng',sprintf('figures\\I_0=0%d-split.png', I_0));

    saveas(fig2, sprintf('figures\\I_0=0%d', I_0),'fig');
    print(fig2,'-dpng',sprintf('figures\\I_0=0%d.png', I_0));
end


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
ylabel('Normal cells (10^{11})', 'fontsize',12)
title(sprintf('Minimum normal cells population = %g', min(Cells_out.data(:,1))), 'fontsize',12)

subplot(4,1,2);
plot(Cells_out.time,Cells_out.data(:,3), 'LineWidth',1)
axis([0 tf 0 1.8])
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Immune cells (10^{11})', 'fontsize',12)
title(sprintf('Initial immune cells population I_0 = 0.%d', I_0), 'fontsize',12)

subplot(4,1,3);
plot(Cells_out.time,Cells_out.data(:,2), 'LineWidth',1)
set(gca,'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Tumor cells (10^{11})', 'fontsize',12)
title(sprintf('Maximum tumor population = %g',  max(Cells_out.data(:,2))), 'fontsize',12)

subplot(4,1,4);
% stairs(t2, v_bb, 'LineWidth',1)
plot(Cells_out.time,Cells_out.data(:,4), 'LineWidth',1)
axis([0 tf 0 1.2])
set(gca, 'FontSize',11)
xlabel('Days', 'fontsize',12)
ylabel('Drug concentration (mg/L)', 'fontsize',12)
title(sprintf('Maximum drug concentration: %g',max(Cells_out.data(:,4))), 'fontsize',12)

% title(sprintf('Total drug : %g mg/L',sum(v_bb)), 'fontsize',12)

%%%%%%%% Figure /w hold
fig2_bangbang = figure();
t2 = linspace(1,tf,tf);
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
ylabel('Cells (10^{11}), Drug (mg/m^2)', 'fontsize',12)
legend('N', 'T', 'I', 'v')

% Save the results
if (save == 1)    
    % Print - All populations & Drug input
    saveas(fig1_bangbang, sprintf('figures\\I_0=0%d-split_bangbang', I_0),'fig');
    print(fig1_bangbang,'-dpng',sprintf('figures\\I_0=0%d-split_bangbang.png', I_0));
    
    % Print - All populations & Drug input
    saveas(fig2_bangbang, sprintf('figures\\I_0=0%d_bangbang', I_0),'fig');
    print(fig2_bangbang,'-dpng',sprintf('figures\\I_0=0%d_bangbang.png', I_0));
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Print-Pulsed!                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% fig_pulsed = figure();
% hold on;
% plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,1), 'LineWidth',1)
% plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,2), 'LineWidth',1)
% plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,3), 'LineWidth',1)
% stairs(pulsed_Cells_out.time, pulsed_Drug_out.data, 'LineWidth',1)
% % plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,4), 'LineWidth',1)
% axis([0 tf 0 1.2])
% set(gca,'FontSize',11)
% title('Traditional Pulsed approach - Cell Populations and Drug input', 'fontsize',12)
% xlabel('Days', 'fontsize',12)
% ylabel('Cells (10^{11}), Drug (mg/m^2)', 'fontsize',12)
% legend('N', 'T', 'I', 'v')
% 
% % Save the results
% if (save == 1)    
%     % Print - All populations & Drug input
%     I_0 = int8(I0*100);
%     saveas(fig_pulsed, sprintf('figures\\I_0=0%d_pulsed', I_0),'fig');
%     print(fig_pulsed,'-dpng',sprintf('figures\\I_0=0%d_pulsed.png', I_0));    
% end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                              Print-DrugFree!                            %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% 
% v_bb_ts = timeseries(0);
% 
% sim('model\\model_depillis_bangbang',simTime);
% 
% 
% fig_drugfree = figure();
% hold on;
% FIX plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,1), 'LineWidth',1)
% FIX plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,2), 'LineWidth',1)
% FIX plot(pulsed_Cells_out.time,pulsed_Cells_out.data(:,3), 'LineWidth',1)
% axis([0 tf 0 1])
% set(gca,'FontSize',11)
% title(sprintf('Drug-free Cell Populations (I_0=0.%d)', I_0), 'fontsize',12)
% xlabel('Days', 'fontsize',12)
% ylabel('Cells (10^{11})', 'fontsize',12)
% legend('N', 'T', 'I')
% 
% % Save the results
% if (save == 1)    
%     % Print - All populations & Drug input    
%     saveas(fig_drugfree, sprintf('figures\\I_0=0%d_drugfree', I_0),'fig');
%     print(fig_drugfree,'-dpng',sprintf('figures\\I_0=0%d_drugfree.png', I_0));    
% end
