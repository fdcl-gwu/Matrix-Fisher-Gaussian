function [ Miu, Sigma, P, U, S, V ] = MFGIMUUpdate( Miu1, Sigma1, P1, U1, S1, V1, MF, z, Sigmaz, bool_prog )

H = [zeros(3,3),eye(3),zeros(3,6)];

if bool_prog
    tau = 0.3;
end

%% attitude
Sigmac1 = Sigma1-P1*(trace(S1)*eye(3)-S1)*P1';
Sigma_R = H*Sigmac1*H'+Sigmaz;

if ~bool_prog
    % select sigma points from prior distribution
    [R1,w1] = pdf_MF_unscented_transform2(U1*S1*V1'+MF);

    % reweight
    for i = 1:7
        vR1 = vee(U1'*R1(:,:,i)*V1*S1-S1*V1'*R1(:,:,i)'*U1);
        Miuc = Miu1+P1*vR1;
        expfR1 = exp(-1/2*(H*Miuc-z)'*Sigma_R^-1*(H*Miuc-z));
        w1(i) = w1(i)*expfR1;
    end
    w1 = w1/sum(w1);
    
    % calculate ER
    ER = sum(R1.*permute(w1,[2,3,1]),3);
else
    lambda_rem = 1;
    expfR1 = zeros(7,1);
    
    [Rp,wp] = pdf_MF_unscented_transform2(U1*S1*V1'+MF);
    
    while lambda_rem > 0
        for i = 1:7
            vR1 = vee(U1'*Rp(:,:,i)*V1*S1-S1*V1'*Rp(:,:,i)'*U1);
            Miuc = Miu1+P1*vR1;
            expfR1(i) = exp(-1/2*lambda_rem*(H*Miuc-z)'*Sigma_R^-1*(H*Miuc-z));
        end
        
        r = min(expfR1)/max(expfR1);
        if r < tau
            k = log(tau)/log(r);
            lambda = lambda_rem*log(tau)/log(r);
            lambda_rem = lambda_rem-lambda;

            expfR1 = expfR1.^k;
        else
            lambda_rem = 0;
        end
        
        wp = wp.*expfR1;
        wp = wp/sum(wp);
        ERp = sum(Rp.*permute(wp,[2,3,1]),3);
        
        if lambda_rem > 0
            [~,Dp,~] = psvd(ERp);
            [Rp,wp] = pdf_MF_unscented_transform2(ERp,[],diag(Dp));
        end
    end
    
    ER = ERp;
end

% update U, S, V
[U,EQ,V] = psvd(ER);
S = diag(pdf_MF_M2S(diag(EQ),diag(S1)));
s = diag(S);

%% linear components
% linear update
Sigmac1 = Sigma1 - P1*(trace(S1)*eye(3)-S1)*P1';
K = Sigmac1*H'*(Sigmaz+H*Sigmac1*H')^-1;

% EQQ
EQQ = zeros(3,3);
for i = 1:3
    for j = setdiff(1:3,i)
        EQQ(3*(i-1)+j,3*(i-1)+j) = EQ(i,i)*s(i)/(s(i)^2-s(j)^2)-EQ(j,j)*s(j)/(s(i)^2-s(j)^2);
        EQQ(3*(i-1)+j,3*(j-1)+i) = EQ(i,i)*s(j)/(s(i)^2-s(j)^2)-EQ(j,j)*s(i)/(s(i)^2-s(j)^2);
    end
end

EQQ(1,1) = 1-EQQ(2,2)-EQQ(3,3);
EQQ(5,5) = 1-EQQ(4,4)-EQQ(6,6);
EQQ(9,9) = 1-EQQ(7,7)-EQQ(8,8);

EQQ(1,5) = EQQ(2,4)+EQ(3,3);
EQQ(1,9) = EQQ(3,7)+EQ(2,2);
EQQ(5,9) = EQQ(6,8)+EQ(1,1);
EQQ(5,1) = EQQ(1,5);
EQQ(9,1) = EQQ(1,9);
EQQ(9,5) = EQQ(5,9);

% EvRvR
EvRvR = zeros(3,3);
EvRvR(1,1) = S(2,2)*EQ(2,2) + S(3,3)*EQ(3,3);
EvRvR(2,2) = S(1,1)*EQ(1,1) + S(3,3)*EQ(3,3);
EvRvR(3,3) = S(1,1)*EQ(1,1) + S(2,2)*EQ(2,2);

% EvR1, EvR1vR1
UT = U1'*U;
VT = V1'*V;
ST = UT'*S1*VT;

EvR1 = UT*vee(EQ*ST'-ST*EQ');

EvRTvRT(1,1) = (ST(3,2)^2)*EQQ(5,5)+(-2*ST(2,3)*ST(3,2))*EQQ(5,9)+(ST(2,3)^2)*EQQ(9,9)...
    +(ST(3,1)^2)*EQQ(2,2)+(ST(2,1)^2)*EQQ(3,3)+(ST(2,2)^2+ST(3,3)^2)*EQQ(6,6)+(-2*ST(2,2)*ST(3,3))*EQQ(6,8);
EvRTvRT(2,2) = (ST(3,1)^2)*EQQ(1,1)+(-2*ST(1,3)*ST(3,1))*EQQ(1,9)+(ST(1,3)^2)*EQQ(9,9)...
    +(ST(3,2)^2)*EQQ(2,2)+(ST(1,1)^2+ST(3,3)^2)*EQQ(3,3)+(-2*ST(1,1)*ST(3,3))*EQQ(3,7)+(ST(1,2)^2)*EQQ(6,6);
EvRTvRT(3,3) = (ST(2,1)^2)*EQQ(1,1)+(-2*ST(1,2)*ST(2,1))*EQQ(1,5)+(ST(1,2)^2)*EQQ(5,5)...
    +(ST(1,1)^2+ST(2,2)^2)*EQQ(2,2)+(-2*ST(1,1)*ST(2,2))*EQQ(2,4)+(ST(2,3)^2)*EQQ(3,3)+(ST(1,3)^2)*EQQ(6,6);
EvRTvRT(1,2) = (-ST(3,1)*ST(3,2))*EQQ(1,5)+(ST(2,3)*ST(3,1))*EQQ(1,9)+(ST(1,3)*ST(3,2))*EQQ(5,9)+(-ST(1,3)*ST(2,3))*EQQ(9,9)...
    +(-ST(3,1)*ST(3,2))*EQQ(2,4)+(-ST(1,1)*ST(2,1))*EQQ(3,3)+(ST(2,1)*ST(3,3))*EQQ(3,7)+(-ST(1,2)*ST(2,2))*EQQ(6,6)+(ST(1,2)*ST(3,3))*EQQ(6,8);
EvRTvRT(1,3) = (ST(2,1)*ST(3,2))*EQQ(1,5)+(-ST(2,1)*ST(2,3))*EQQ(1,9)+(-ST(1,2)*ST(3,2))*EQQ(5,5)+(ST(1,2)*ST(2,3))*EQQ(5,9)...
    +(-ST(1,1)*ST(3,1))*EQQ(2,2)+(ST(2,2)*ST(3,1))*EQQ(2,4)+(-ST(2,1)*ST(2,3))*EQQ(3,7)+(-ST(1,3)*ST(3,3))*EQQ(6,6)+(ST(1,3)*ST(2,2))*EQQ(6,8);
EvRTvRT(2,3) = (-ST(2,1)*ST(3,1))*EQQ(1,1)+(ST(1,2)*ST(3,1))*EQQ(1,5)+(ST(1,3)*ST(2,1))*EQQ(1,9)+(-ST(1,2)*ST(1,3))*EQQ(5,9)...
    +(-ST(2,2)*ST(3,2))*EQQ(2,2)+(ST(1,1)*ST(3,2))*EQQ(2,4)+(-ST(2,3)*ST(3,3))*EQQ(3,3)+(ST(1,1)*ST(2,3))*EQQ(3,7)+(-ST(1,2)*ST(1,3))*EQQ(6,8);
EvRTvRT(2,1) = EvRTvRT(1,2);
EvRTvRT(3,1) = EvRTvRT(1,3);
EvRTvRT(3,2) = EvRTvRT(2,3);

EvR1vR1 = UT*EvRTvRT*UT';

% EvR1vR
EvRTvR(1,1) = (ST(2,2)*s(2)+ST(3,3)*s(3))*EQQ(6,6) +(-ST(2,2)*s(3)-ST(3,3)*s(2))*EQQ(6,8);
EvRTvR(1,2) = (-ST(2,1)*s(1))*EQQ(3,3) +(ST(2,1)*s(3))*EQQ(3,7);
EvRTvR(1,3) = (-ST(3,1)*s(1))*EQQ(2,2) +(ST(3,1)*s(2))*EQQ(2,4);
EvRTvR(2,1) = (-ST(1,2)*s(2))*EQQ(6,6) +(ST(1,2)*s(3))*EQQ(6,8);
EvRTvR(2,2) = (ST(1,1)*s(1)+ST(3,3)*s(3))*EQQ(3,3) +(-ST(1,1)*s(3)-ST(3,3)*s(1))*EQQ(3,7);
EvRTvR(2,3) = (-ST(3,2)*s(2))*EQQ(2,2) +(ST(3,2)*s(1))*EQQ(2,4);
EvRTvR(3,1) = (-ST(1,3)*s(3))*EQQ(6,6) +(ST(1,3)*s(2))*EQQ(6,8);
EvRTvR(3,2) = (-ST(2,3)*s(3))*EQQ(3,3) +(ST(2,3)*s(1))*EQQ(3,7);
EvRTvR(3,3) = (ST(1,1)*s(1)+ST(2,2)*s(2))*EQQ(2,2) +(-ST(1,1)*s(2)-ST(2,2)*s(1))*EQQ(2,4);

EvR1vR = UT*EvRTvR;

% Ex, Exx, ExvR
IKH = eye(12)-K*H;

Ex = IKH*(Miu1+P1*EvR1) + K*z;
Exx = IKH*Sigmac1 + IKH*(Miu1*Miu1'+Miu1*EvR1'*P1'+P1*EvR1*Miu1'+P1*EvR1vR1*P1')*IKH' +...
    IKH*(Miu1+P1*EvR1)*z'*K' + K*z*(Miu1'+EvR1'*P1')*IKH' + K*(z*z')*K';
ExvR = IKH*P1*EvR1vR;

% covxx, covxvR, covvRvR
covxx = Exx-Ex*Ex';
covxvR = ExvR;
covvRvR = EvRvR;

%% estimate
P = covxvR*covvRvR^-1;
Miu = Ex;
Sigma = covxx-P*covxvR'+P*(trace(S)*eye(3)-S)*P';

end

