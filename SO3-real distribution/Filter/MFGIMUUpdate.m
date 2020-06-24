function [ Miu, Sigma, P, U, S, V ] = MFGIMUUpdate( Miu1, Sigma1, P1, U1, S1, V1, MF, z, zSigma)

%% attitude
% select sigma points from prior distribution
[R1,w1] = MFGetSigmaPoints(U1*S1*V1'+MF);

% calculate E[R]
Sigmac1 = Sigma1-P1*(trace(S1)*eye(3)-S1)*P1';

ER = zeros(3);
expfR1 = zeros(1,7);
vR1 = zeros(3,7);
for i = 1:7
    vR1(:,i) = vee(U1'*R1(:,:,i)*V1*S1-S1*V1'*R1(:,:,i)'*U1);
    expfR1(i) = exp(-1/2*(Miu1+P1*vR1(:,i)-z)'*(Sigmac1+zSigma)^-1*(Miu1+P1*vR1(:,i)-z));
    ER = ER + w1(i)*R1(:,:,i)*expfR1(i);
end

% update U, S, V
[U,D,V] = psvd(ER);
S = diag(pdf_MF_M2S(diag(D),diag(S1)));

%% linear components
% EvR, EvRvR
EvR = zeros(3,1);
EvRvR = zeros(3);
vR = zeros(3,7);
for i = 1:7
    vR(:,i) = vee(U'*R1(:,:,i)*V*S-S*V'*R1(:,:,i)'*U);
    EvR = EvR + w(i)*vR(:,i)*expfR1(i);
    EvRvR = EvRvR + w(i)*(vR(:,i)*vR(:,i)')*expfR1(i);
end

% EvR1, EvR1vR, EvR1vR1
EvR1 = zeros(3,1);
EvR1vR = zeros(3,3);
EvR1vR1 = zeros(3,3);
for i = 1:7
    EvR1 = EvR1 + w1(i)*vR1(:,i)*expfR1(i);
    EvR1vR = EvR1vR + w1(i)*vR1(:,i)*vR(:,i)'*expfR1(i);
    EvR1vR1 = EvR1vR1 + w1(i)*vR1(:,i)*vR1(:,i)'*expfR1(i);
end

% Ex, ExvR, Exx
MiuT = zSigma*(Sigmac1+zSigma)^-1*Miu1 + Sigmac1*(Sigmac1+zSigma)^-1*z;
PT = zSigma*(Sigmac1+zSigma)^-1*P1;

Ex = PT*EvR1 + MiuT;
ExvR = PT*EvR1vR + MiuT*EvR;
Exx = (Sigmac1^-1+zSigma^-1)^-1 + MiuT*MiuT' + MiuT*EvR1'*PT' + PT*EvR1*MiuT' + PT*EvR1vR1*PT';

% covxx, covxvR, covvRvR
covxx = Exx-Ex*Ex';
covxvR = ExvR-Ex*EvR';
covvRvR = EvRvR - EvR*EvR';

% estimate
P = covxvR*covvRvR^-1;
Miu = Ex-P*EvR;
Sigma = covxx-P*covxvR'+P*(trace(S)*eye(3)-S)*P';

end

