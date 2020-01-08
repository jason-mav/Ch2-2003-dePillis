close all;
clear all;

parameters;

step = 0.05;


%% Figure 1
% axis
% [I, N] = meshgrid(-10:step:10,-10:step:10);    
% [I, T] = meshgrid(-10:step:10,-10:step:10);

%% NN
[I, T] = meshgrid(0:step:9,0:step:4);  

N = 1/b2 - (c4/(r2*b2)).*T;

surface(N,T,I)
axis([0 1.5 0 3.5 0 3])
xlabel('N');
ylabel('T');
zlabel('I');

%% NI
[T, N] = meshgrid(0:step:3.7,0:step:1.5);

I = s.*(alpha+T)./(c1.*T.*(alpha+T) + d1.*(alpha+T) - ro.*T);

surface(N,T,I)

%     hold on;

% axis([0 1.5 0 3.7 0 9]);
xlabel('N');
ylabel('T');
zlabel('I');

% %% Figure 1 - N2 & N3
% % axis
% [I, N] = meshgrid(0:step:9,0:step:1.5);  


%% NT
[I, N] = meshgrid(0:step:9,0:step:1.5);

T = 1/b1 - (c2/(r1*b1)).*I - (c3/(r1*b1)).*N;

surface(N,T,I)

% axis([0 1.5 0 3.7 0 9])
xlabel('Normal');
ylabel('Tumor');
zlabel('Immune');
hold on;


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    