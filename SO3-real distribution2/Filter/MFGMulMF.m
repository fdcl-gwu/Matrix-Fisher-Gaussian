function [ Miu, Sigma, P, U, S, V ] = MFGMulMF( Miu1, Sigma1, P1, U1, S1, V1, FM )

%% Matrix Fisher part
% parameters
[U,S,V] = psvd(U1*S1*V1'+FM);

% moments
EQ = pdf_MF_moment(diag(S));
EQQ = pdf_MF_moment2(diag(S));

%% fR1
% first moment
EQ1 = U1'*U*diag(EQ)*V'*V1;
EfR1 = vee(S1*EQ1-EQ1'*S1);

% second moment
dU = U1'*U;
dV = V1'*V;
dS = dU'*S1*dV;
EfRfR1 = dV*EfRfRX(EQQ,dS)*dV';

%% fR2
EfR2 = [0;0;0];
EfRfR2 = EfRfRX(EQQ,S);

%% fR1 fR2
EfRfR12 = dV*EfRfRXS(EQQ,dS,diag(S));

%% estimation
SigmaMInv1 = diag([S1(2,2)+S1(3,3),S1(1,1)+S1(3,3),S1(1,1)+S1(2,2)]);
Sigmac1 = Sigma1-P1*SigmaMInv1*P1';
SigmaMInv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);

Ex = Miu1+P1*EfR1;
Exx = Miu1*Miu1'+Miu1*EfR1'*P1'+P1*EfR1*Miu1'+P1*EfRfR1*P1'+Sigmac1;
ExfR2 = Miu1*EfR2'+P1*EfRfR12;

covxx = Exx-Ex*Ex';
covxfR2 = ExfR2-Ex*EfR2';
covfRfR2 = EfRfR2-EfR2*EfR2';

P = covxfR2*covfRfR2^-1;
Miu = Ex-P*EfR2;
Sigma = covxx-P*covxfR2'+P*SigmaMInv*P';

end


function [ E ] = EfRfRX( EQQ, ST )

E(1,1) = (ST(2,3)^2)*EQQ(5,5)+(-2*ST(2,3)*ST(3,2))*EQQ(5,9)...
    +(ST(3,2)^2)*EQQ(9,9)+(ST(1,3)^2)*EQQ(2,2)...
    +(ST(1,2)^2)*EQQ(3,3)+(ST(2,2)^2+ST(3,3)^2)*EQQ(6,6)+(-2*ST(2,2)*ST(3,3))*EQQ(6,8);
E(2,2) = (ST(1,3)^2)*EQQ(1,1)+(-2*ST(1,3)*ST(3,1))*EQQ(1,9)...
    +(ST(3,1)^2)*EQQ(9,9)+(ST(2,3)^2)*EQQ(2,2)...
    +(ST(1,1)^2+ST(3,3)^2)*EQQ(3,3)+(-2*ST(1,1)*ST(3,3))*EQQ(3,7)+(ST(2,1)^2)*EQQ(6,6);
E(3,3) = (ST(1,2)^2)*EQQ(1,1)+(-2*ST(1,2)*ST(2,1))*EQQ(1,5)...
    +(ST(2,1)^2)*EQQ(5,5)+(ST(1,1)^2+ST(2,2)^2)*EQQ(2,2)...
    +(-2*ST(1,1)*ST(2,2))*EQQ(2,4)+(ST(3,2)^2)*EQQ(3,3)+(ST(3,1)^2)*EQQ(6,6);
E(1,2) = (-ST(1,3)*ST(2,3))*EQQ(1,5)+(ST(1,3)*ST(3,2))*EQQ(1,9)...
    +(ST(2,3)*ST(3,1))*EQQ(5,9)+(-ST(3,1)*ST(3,2))*EQQ(9,9)...
    +(-ST(1,3)*ST(2,3))*EQQ(2,4)+(-ST(1,1)*ST(1,2))*EQQ(3,3)...
    +(ST(1,2)*ST(3,3))*EQQ(3,7)+(-ST(2,1)*ST(2,2))*EQQ(6,6)+(ST(2,1)*ST(3,3))*EQQ(6,8);
E(1,3) = (ST(1,2)*ST(2,3))*EQQ(1,5)+(-ST(1,2)*ST(3,2))*EQQ(1,9)...
    +(-ST(2,1)*ST(2,3))*EQQ(5,5)+(ST(2,1)*ST(3,2))*EQQ(5,9)...
    +(-ST(1,1)*ST(1,3))*EQQ(2,2)+(ST(1,3)*ST(2,2))*EQQ(2,4)...
    +(-ST(1,2)*ST(3,2))*EQQ(3,7)+(-ST(3,1)*ST(3,3))*EQQ(6,6)+(ST(2,2)*ST(3,1))*EQQ(6,8);
E(2,3) = (-ST(1,2)*ST(1,3))*EQQ(1,1)+(ST(1,3)*ST(2,1))*EQQ(1,5)...
    +(ST(1,2)*ST(3,1))*EQQ(1,9)+(-ST(2,1)*ST(3,1))*EQQ(5,9)...
    +(-ST(2,2)*ST(2,3))*EQQ(2,2)+(ST(1,1)*ST(2,3))*EQQ(2,4)...
    +(-ST(3,2)*ST(3,3))*EQQ(3,3)+(ST(1,1)*ST(3,2))*EQQ(3,7)+(-ST(2,1)*ST(3,1))*EQQ(6,8);
E(2,1) = E(1,2);
E(3,1) = E(1,3);
E(3,2) = E(2,3);

end


function [ E ] = EfRfRXS( EQQ, ST, s )

E(1,1) = (ST(2,2)*s(2)+ST(3,3)*s(3))*EQQ(6,6) + (-ST(2,2)*s(3)-ST(3,3)*s(2))*EQQ(6,8);
E(1,2) = (-ST(1,2)*s(1))*EQQ(3,3) + (ST(1,2)*s(3))*EQQ(3,7);
E(1,3) = (-ST(1,3)*s(1))*EQQ(2,2) + (ST(1,3)*s(2))*EQQ(2,4);
E(2,1) = (-ST(2,1)*s(2))*EQQ(6,6) + (ST(2,1)*s(3))*EQQ(6,8);
E(2,2) = (ST(1,1)*s(1)+ST(3,3)*s(3))*EQQ(3,3) + (-ST(1,1)*s(3)-ST(3,3)*s(1))*EQQ(3,7);
E(2,3) = (-ST(2,3)*s(2))*EQQ(2,2) + (ST(2,3)*s(1))*EQQ(2,4);
E(3,1) = (-ST(3,1)*s(3))*EQQ(6,6) + (ST(3,1)*s(2))*EQQ(6,8);
E(3,2) = (-ST(3,2)*s(3))*EQQ(3,3) + (ST(3,2)*s(1))*EQQ(3,7);
E(3,3) = (ST(1,1)*s(1)+ST(2,2)*s(2))*EQQ(2,2) + (-ST(1,1)*s(2)-ST(2,2)*s(1))*EQQ(2,4);

end

