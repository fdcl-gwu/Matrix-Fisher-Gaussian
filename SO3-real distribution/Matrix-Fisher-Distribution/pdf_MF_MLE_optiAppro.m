function [ F ] = pdf_MF_MLE_optiAppro( ER, S0, U0, V0 )

% default initialization
if ~exist('S0','var')
    S0 = diag([1,1,1]);
    U0 = eye(3);
    V0 = eye(3);
end

% step size
kS = 1;
kR = 0.01;

% first step
c0 = pdf_MF_normal(diag(S0));
fOld = log(c0)-trace(V0*S0*U0'*ER);
[dS,dU,dV,c] = diffLogLike(ER,S0,U0,V0);
S = S0-kS*dS;
U = U0*expRM(-kR*dU);
V = V0*expRM(-kR*dV);
f = log(c)-trace(V*S*U'*ER);

% iteration
while abs(f-fOld)>1e-8
    while abs(f-fOld)>1e-6 && sqrt(sum(dV.^2))>1e-3
        while abs(f-fOld)>1e-6 && sqrt(sum(dU.^2))>1e-3
            fOld = f;
            [dS,dU,dV,c] = diffLogLike(ER,S,U,V);
            U = U*expRM(-kR*dU);
            f = log(c)-trace(V*S*U'*ER)
        end
        
        fOld = f;
        [dS,dU,dV,c] = diffLogLike(ER,S,U,V);
        V = V*expRM(-kR*dV);
        f = log(c)-trace(V*S*U'*ER)
    end
    
    fOld = f;
    [dS,dU,dV,c] = diffLogLike(ER,S,U,V);
    S = S-kS*dS;
    F = U*S*V';
    f = log(c)-trace(V*S*U'*ER)
end

end


function [ dS, dU, dV, c ] = diffLogLike( R, S, U, V )

% with respect to S
c = pdf_MF_normal(diag(S));
dc = pdf_MF_normal_deriv(diag(S));
dS1 = diag(dc/c);
Q = U'*R*V;
dS2 = -diag([Q(1,1),Q(2,2),Q(3,3)]);
dS = dS1+dS2;

% with respect to U
QS = Q*S;
dU(1,1) = QS(2,3)-QS(3,2);
dU(2,1) = -QS(1,3)+QS(3,1);
dU(3,1) = QS(1,2)-QS(2,1);

% with respect to V
SQ = S*Q;
dV(1,1) = -SQ(2,3)+SQ(3,2);
dV(2,1) = +SQ(1,3)-SQ(3,1);
dV(3,1) = -SQ(1,2)+SQ(2,1);

end

