syms r1 r2 b1 b2 d1 d2 c1 c2 c3 c4 N T I a1 a2 a3 s ro alpha v u x1 x2 x3 x4;

N_dot = r2*N*(1 -b2*N) -c4*T*N -a3*(1-exp(-u))*N;
N_dot_int = int(N_dot)

T_dot = r1*T*(1 -b1*T) -c2*I*T -c3*T*N -a2*(1-exp(-u))*T;
T_dot_int = int(T_dot)

I_dot = s +ro*I*T/(alpha +T) -c1*I*T -d1*I -a1*(1-exp(-u))*I;
I_dot_int = int(I_dot)

u_dot = v - d2*u;
u_dot_int = int(u_dot)