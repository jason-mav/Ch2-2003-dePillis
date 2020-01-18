% function [I, N2, N3] = Nullspaces(s, alpha, ro, c1, c2, c3, c4, d1, b1, r1, r2, I, T, N)
%     c1 = 0.558;   %0.558 %0.6;
%     c2 = 0.4875;   %39/80;
%     c3 = 1.0125;     %81/80;
%     r1 = 1.5;   %c2+c3
%     ro = 1; %1;

close all;
clear all;
parameters;

step = 0.05;


%% Figure 1 - N1
% axis
[I, N] = meshgrid(0:step:5,2:-step:0);    
[I, T] = meshgrid(0:step:5,0:step:2);

% N1
[T, N] = meshgrid(0:step:2,2:-step:0);

I = s.*(alpha+T)./(c1.*T.*(alpha+T) + d1.*(alpha+T) - ro.*T);

surface(N,I,T)
%     hold on;

axis([0 2 -5 5 -1 2]);
xlabel('N');
ylabel('I');
zlabel('T');

%% Figure 1 - N2 & N3
% axis
[I, N] = meshgrid(0:step:2,2:-step:0);  

% N2
figure;
T = 1/b1 - (c2/r1*b1).*I - (c3/r1*b1).*N;

surface(N,I,T)

axis([-1 2 0 2 -1 2])
xlabel('N');
ylabel('I');
zlabel('T');
hold on;

% N3
[I, T] = meshgrid(0:step:2,0:step:2);  

N = 1 - (c4/r2).*T;

surface(N,I,T)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    