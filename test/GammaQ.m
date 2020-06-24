syms ST11 ST12 ST13 ST21 ST22 ST23 ST31 ST32 ST33
syms Q11 Q12 Q13 Q21 Q22 Q23 Q31 Q32 Q33

GammaQ(1,1) = ST22*Q33+ST33*Q22-ST23*Q32-ST32*Q23;
GammaQ(1,2) = ST23*Q31+ST31*Q23-ST21*Q33-ST33*Q21;
GammaQ(1,3) = ST21*Q32+ST32*Q21-ST22*Q31-ST31*Q22;
GammaQ(2,1) = ST13*Q32+ST32*Q13-ST12*Q33-ST33*Q12;
GammaQ(2,2) = ST11*Q33+ST33*Q11-ST13*Q31-ST31*Q13;
GammaQ(2,3) = ST12*Q31+ST31*Q12-ST11*Q32-ST32*Q11;
GammaQ(3,1) = ST12*Q23+ST23*Q12-ST13*Q22-ST22*Q13;
GammaQ(3,2) = ST13*Q21+ST21*Q13-ST11*Q23-ST23*Q11;
GammaQ(3,3) = ST11*Q22+ST22*Q11-ST12*Q21-ST21*Q12;
