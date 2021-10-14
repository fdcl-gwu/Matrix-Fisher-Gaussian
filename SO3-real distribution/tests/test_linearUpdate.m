close all;
clear;

addpath('../Matrix-Fisher-Distribution');
addpath('../../rotation3d');

% prior density
S = diag([10,10,10]);
U = eye(3);
V = eye(3);
Miu = 0;
Sigma = 1;
P = [0.7/sqrt(20),0,0];

% measurement
Mium = 5;
Sigmam = 5;

% marginal likelihood
Sigmac = Sigma-P*(trace(S)*eye(3)-S)*P';
Sigma_post = Sigmam+Sigmac;

% Monte Carlo method
Ns = 10000;
R = pdf_MF_sampling(S,Ns);
w = 1/Ns*ones(Ns,1);

for i = 1:Ns
    vR = vee(R(:,:,i)*S-S*R(:,:,i)');
    l = exp(-0.5*(Miu+P*vR-Mium)^2/Sigma_post);
    w(i) = w(i)*l;
end
w = w/sum(w);

ER = sum(R.*permute(w,[2,3,1]),3);
[Un_MC,Dn_MC,Vn_MC] = psvd(ER);
Sn_MC = diag(pdf_MF_M2S(diag(Dn_MC)));

% sigma points method
[R,w] = pdf_MF_unscented_transform2(S);

for i = 1:7
    vR = vee(R(:,:,i)*S-S*R(:,:,i)');
    l = exp(-0.5*(Miu+P*vR-Mium)^2/Sigma_post);
    w(i) = w(i)*l;
end
w = w/sum(w);

ER = sum(R.*permute(w,[2,3,1]),3);
[Un_UT,Dn_UT,Vn_UT] = psvd(ER);
Sn_UT = diag(pdf_MF_M2S(diag(Dn_UT)));

% progressive update
[Rp,wp] = pdf_MF_unscented_transform2(S);

tau = 0.3;
l = zeros(7,1);
lambda_rem = 1;

while lambda_rem ~= 0
    for i = 1:7
        Q = U'*Rp(:,:,i)*V;
        vR = vee(Q*S-S*Q');
        l(i) = exp(-0.5*lambda_rem*(Miu+P*vR-Mium)^2/Sigma_post);
    end
    
    r = min(l)/max(l);
    if r < tau
        k = log(tau)/log(r);
        lambda = lambda_rem*log(tau)/log(r);
        lambda_rem = lambda_rem-lambda;
        
        l = l.^k;
    else
        lambda_rem = 0;
    end
    
    wp = wp.*l;
    wp = wp/sum(wp);
    ERp = sum(Rp.*permute(wp,[2,3,1]),3);
    [Up,Dp,Vp] = psvd(ERp);
    [Rp,wp] = pdf_MF_unscented_transform2(ERp,[],diag(Dp));
end

ERp = sum(Rp.*permute(wp,[2,3,1]),3);
[Un_p,Dn_p,Vn_p] = psvd(ERp);
Sn_p = diag(pdf_MF_M2S(diag(Dn_p)));

rmpath('../Matrix-Fisher-Distribution');
rmpath('../../rotation3d');
