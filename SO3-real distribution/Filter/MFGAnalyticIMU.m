function [ R, x, MFG, stepT ] = MFGAnalyticIMU( gyro, acce, RMea, pMea, defQS, parameters )

N = size(gyro,2);
dt = parameters.dt;
sf_GPS = parameters.sf_GPS;

% noise parameters
randomWalk = parameters.randomWalk;
biasInstability = parameters.biasInstability;
acceRandomWalk = parameters.acceRandomWalk;
acceBiasInstability = parameters.acceBiasInstability;
if parameters.GaussMea
    SM = Gau2MF(parameters.rotMeaNoise);
else
    SM = parameters.rotMeaNoise;
end
posMeaNoise = parameters.posMeaNoise;

useGrav = parameters.useGrav;
if useGrav
    gravMeaNoise = Gau2VM(parameters.gravMeaNoise);
end

H.Hgu = randomWalk*eye(3);
H.Hgv = biasInstability*eye(3);
H.Hau = acceRandomWalk*eye(3);
H.Hav = acceBiasInstability*eye(3);

hasRMea = ~isempty(RMea);
bool_prog = parameters.bool_prog;

noGPS = parameters.noGPS;

% options for lsqlin
options = optimoptions('lsqlin','display','off','algorithm','active-set');

% initialize distribution
Miu = parameters.xInit;
Sigma = parameters.initXNoise;
P = zeros(12,3);
if parameters.GaussMea
    U = parameters.RInit;
    V = eye(3);
    S = Gau2MF(parameters.initRNoise);
else
    [U,S,V] = psvd(parameters.initRNoise);
end
S(2,2) = S(2,2)+1e-5;
S(3,3) = S(3,3)+2e-5;

% data containers
MFG.Miu = zeros(12,N); MFG.Miu(:,1) = Miu;
MFG.Sigma = zeros(12,12,N); MFG.Sigma(:,:,1) = Sigma;
MFG.P = zeros(12,3,N); MFG.P(:,:,1) = P;
MFG.U = zeros(3,3,N); MFG.U(:,:,1) = U;
MFG.V = zeros(3,3,N); MFG.V(:,:,1) = V;
MFG.S = zeros(3,N); MFG.S(:,1) = diag(S);
R = zeros(3,3,N); R(:,:,1) = U*V';
x = zeros(12,N); x(:,1) = Miu;
stepT = zeros(N-1,1);

% filter iteration
for n = 2:N
    tic;
    
    % uncertainty propagation
    omega = (gyro(:,n-1)+gyro(:,n))/2;
    a = (acce(:,n-1)+acce(:,n))/2;
    [Miu,Sigma,P,U,S,V] = MFGIMUProp(omega,a,Miu,Sigma,P,U,S,V,H,defQS,dt,options);
    
    % update
    if rem(n,1/dt/sf_GPS)==0 && (isempty(noGPS) || n<noGPS(1) || n>noGPS(2))
        if hasRMea
            if useGrav
                grav = acce(:,n)+Miu(10:12);
                FMea = RMea(:,:,n)*SM + gravMeaNoise*[0;0;1]*grav'/sqrt(sum(grav.^2));
            else
                FMea = RMea(:,:,n)*SM;
            end
        else
            if useGrav
                grav = acce(:,n)+Miu(10:12);
                FMea = gravMeaNoise*[0;0;1]*grav'/sqrt(sum(grav.^2));
            else
                FMea = zeros(3,3);
            end
        end
        
        [Miu,Sigma,P,U,S,V] = MFGIMUUpdate(Miu,Sigma,P,U,S,V,FMea,pMea(:,n),posMeaNoise,defQS,bool_prog,options);
    else
        if useGrav
            grav = acce(:,n)+Miu(10:12);
            FMea = gravMeaNoise*[0;0;1]*grav'/sqrt(sum(grav.^2));
            [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,FMea,defQS);
        end
    end
    
    % record results
    MFG.Miu(:,n) = Miu;
    MFG.Sigma(:,:,n) = Sigma;
    MFG.P(:,:,n) = P;
    MFG.U(:,:,n) = U;
    MFG.V(:,:,n) = V;
    MFG.S(:,n) = diag(S);
    R(:,:,n) = U*V';
    x(:,n) = Miu;
    
    stepT(n-1) = toc;
end

end


function [ S ] = Gau2MF( Sigma )

N = 100000;
v = mvnrnd([0;0;0],Sigma,N);

R = expRot(v);
ER = mean(R,3);
[~,D,~] = psvd(ER);

S = diag(pdf_MF_M2S(diag(D)));

end


function [ kappa ] = Gau2VM( Sigma )

N = 100000;
v = sqrtm(Sigma)*randn(3,N)+[0;0;9.8];
v = v./sqrt(sum(v.^2));

rho = sqrt(sum(mean(v,2).^2));

options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-15,'Display','off');
kappa = fsolve(@(k) coth(k)-1/k-rho,1,options);

end

