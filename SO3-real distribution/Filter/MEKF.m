function [ R, x, G, stepT ] = MEKF( gyro, mea, parameters )

N = size(gyro,2);
dt = parameters.dt;

% noise parameters
randomWalk = parameters.randomWalk;
biasInstability = parameters.biasInstability;
if parameters.meaIsVec
    vecMeaNoise = parameters.meaNoise;
else
    if parameters.GaussMea
        rotMeaNoise = parameters.rotMeaNoise;
    else
        rotMeaNoise = MF2Gau(parameters.rotMeaNoise);
    end
end

% measurement
if parameters.meaIsVec
    vMea = mea{1};
    vRef = mea{2};
end

% initialize distribution
if parameters.GaussMea
    Sigma = [parameters.initRNoise,zeros(3)
        zeros(3),parameters.initXNoise];
else
    Sigma = [MF2Gau(parameters.initRNoise),zeros(3)
        zeros(3),parameters.initXNoise];
end

% data containers
G.Sigma = zeros(6,6,N); G.Sigma(:,:,1) = Sigma;
R = zeros(3,3,N); R(:,:,1) = parameters.RInit;
x = zeros(3,N); x(:,1) = parameters.xInit;
stepT = zeros(N-1,1);

% filter iteration
for n = 2:N
    tic;
    
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2-x(:,n-1);
    Rp = R(:,:,n-1)*expRot(av*dt);
    F = [expRot(av*dt)',-eye(3)*dt;zeros(3),eye(3)];
    Sigma = F*Sigma*F'+[eye(3)*randomWalk^2*dt,zeros(3);
        zeros(3),eye(3)*biasInstability^2*dt];
    
    % update
    if rem(n,5)==0
        if parameters.meaIsVec
            if size(vRef,1)==3
                H = [hat(Rp'*vRef(:,n)),zeros(3)];
                K = Sigma*H'*(H*Sigma*H'+eye(3)/vecMeaNoise)^-1;
                dx = K*(vMea(:,n)-Rp'*vRef(:,n));
            else
                H = [hat(Rp'*vRef(1:3,n)),zeros(3);hat(Rp'*vRef(4:6,n)),zeros(3)];
                K = Sigma*H'*(H*Sigma*H'+eye(6)/vecMeaNoise)^-1;
                dx = K*[vMea(1:3,n)-Rp'*vRef(1:3,n);vMea(4:6,n)-Rp'*vRef(4:6,n)];
            end
            Sigma = (eye(6)-K*H)*Sigma;
            
            R(:,:,n) = Rp*expRot(dx(1:3));
            x(:,n) = x(:,n-1)+dx(4:6);
        else
            H = [eye(3),zeros(3)];
            K = Sigma*H'*(H*Sigma*H'+rotMeaNoise)^-1;
            dx = K*logRot(Rp'*mea(:,:,n),'v');
            Sigma = (eye(6)-K*H)*Sigma;

            R(:,:,n) = Rp*expRot(dx(1:3));
            x(:,n) = x(:,n-1)+dx(4:6);
        end
    else
        R(:,:,n) = Rp;
        x(:,n) = x(:,n-1);
    end
    
    % record covariance
    G.Sigma(:,:,n) = Sigma;
    
    stepT(n-1) = toc;
end

end


function [ Sigma ] = MF2Gau( S )

N = 100000;
R = pdf_MF_sampling(S,N);

v = logRot(R,'v');
Sigma = cov(v');

end

