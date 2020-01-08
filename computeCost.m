function [vals,derivs] = myCostFunc(params)
% Extract the current design variable values from the parameter object, params.
x = params.Value;



I_dot = s + ro*I*T/(alpha+T) - c1*I*T - d1*I - a1*(1 -exp(-u))*I;
T_dot = r1*T*(1-b1*T) - c2*I*T - c3*T*N - a2*(1 -exp(-u))*T;
N_dot = r2*N*(1-b2*N) - c4*T*N - a3*(1 -exp(-u))*N;
u_dot = v - d2*u;


% Cost function weights
w1 = 1500;
w2 = 150;
w3 = 1000;

y1 = T_data.data(end); % Tumor size at final time
y2 = sum(T_data.data); % total sum of Tumor size
y3 = max(T_data.data); % maximum Tumor size


% Compute the requirements (objective and constraint violations) and 
% assign them to vals, the output of the cost function. 
vals.F = x.^2;
vals.Cleq = x.^2-4*x+1;

% Compute the cost and constraint derivatives.
derivs.F = 2*x;
derivs.Cleq = 2*x-4;

end