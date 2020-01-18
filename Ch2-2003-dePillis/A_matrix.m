function A = A_matrix(a1, a2, a3, b1, b2, c1, c2, c3, c4, d1, d2, r1, r2, ro, s, alpha, x1, x2, x3, x4)

    A11 =  -r2*(1+b2*x1); 
    A12 =  -c4*(x1 + 1/b2);
    A13 =  0; 
    A14 =  -a3*(1/b2 + x1);

    A21 =  -c3*x2;
    A22 =  r1*(1-b1*x2)-(c2*s/d1 + c3/b2);
    A23 =  -c2*x2;
    A24 =  -a2*x2;
        
    A31 =  0;
    A32 =  ro*(x3 + s/d1)/(alpha + x2) -c1*(x3+s/d1) -x4;
    A33 =  -d1; 
    A34 =  -a1*(x3 + s/d1)+x2;
    
    A41 =  0;
    A42 =  0;
    A43 =  0;
    A44 =  -d2;
    
    A = [A11 A12 A13 A14;
         A21 A22 A23 A24;
         A31 A32 A33 A34;
         A41 A42 A43 A44];
end

% syms a1 a2 a3 b1 b2 c1 c2 c3 c4 d1 d2 r1 r2 ro s alpha x1 x2 x3 x4

% x1 = N-1/b2;
% x2 = T;
% x3 = I-s/d1;
% x4 = M; 