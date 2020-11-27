function [ R, x, Sigma, stepT ] = MEKF( gyro, Mea, sf, parameters )

N = size(gyro,2);
dt = 1/sf;

%% settings
% angular velocity noise parameters
if exist('parameters','var') && isfield(parameters,'omegaNoise')
    randomWalk = parameters.omegaNoise.randomWalk;
    biasInstability = parameters.omegaNoise.biasInstability;
else
    randomWalk = 10*pi/180;
    biasInstability = 500/3600*pi/180;
end

% measurement noise parameters
if exist('parameters','var') && isfield(parameters,'meaNoise')
    meaNoise = parameters.meaNoise;
else
    meaNoise = 0.2^2*eye(3);
end

% other settings
if exist('parameters','var') && isfield(parameters,'setting')
    omegaLocal = parameters.setting.omegaLocal;
    GaussMea = parameters.setting.GaussMea;
    meaIsVec = parameters.setting.meaIsVec;
else
    omegaLocal = true;
    GaussMea = true;
    meaIsVec = false;
end

if meaIsVec
    if exist('parameters','var') && isfield(parameters,'setting')
        vecRefInertial = parameters.setting.vecRefInertial;
        nVecRef = parameters.setting.nVecRef;
        vRef = parameters.setting.vRef;
    else
        vecRefInertial = true;
        nVecRef = 1;
        vRef = [0;0;1];
    end
else
    if exist('parameters','var') && isfield(parameters,'setting')
        attMeaLocal = parameters.setting.attMeaLocal;
    else
        attMeaLocal = true;
    end
end

% convert noise distributions
if meaIsVec
    if ~GaussMea
        for nv = 1:nVecRef
            meaNoise(nv) = VM2Gau(meaNoise(nv));
        end
    end
else
    if ~GaussMea
        meaNoise = MF2Gau(meaNoise);
    end
end

%% initialization
% data containers
R = zeros(3,3,N);
x = zeros(3,N);
Sigma = zeros(6,6,N);
stepT = zeros(N-1,1);

% initialize
if exist('parameters','var') && isfield(parameters,'initValue')
    R(:,:,1) = parameters.initValue.U*parameters.initValue.V';
    x(:,1) = parameters.initValue.Miu;
    if GaussMea
        Sigma(:,:,1) = [parameters.initValue.RNoise,zeros(3);
            zeros(3),parameters.initValue.xNoise];
    else
        Sigma(:,:,1) = [MF2Gau(parameters.initValue.RNoise),zeros(3);
            zeros(3),parameters.initValue.xNoise];
    end
else
    R(:,:,1) = eye(3);
    x(:,1) = [0;0;0];
    Sigma(:,:,1) = [0.2^2*eye(3),zeros(3);zeros(3),0.05^2*eye(3)];
end

%% filter iteration
for n = 2:N
    tic;
    
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2-x(:,n-1);
    if omegaLocal
        R(:,:,n) = R(:,:,n-1)*expRot(av*dt);
        F = [expRot(av*dt)',-eye(3)*dt;zeros(3),eye(3)];
        Sigma(:,:,n) = F*Sigma(:,:,n-1)*F'+[eye(3)*randomWalk^2*dt,zeros(3);
            zeros(3),eye(3)*biasInstability^2*dt];
    else
        R(:,:,n) = expRot(av*dt)*R(:,:,n-1);
        F = [eye(3),-R(:,:,n)'*dt;zeros(3),eye(3)];
        Sigma(:,:,n) = F*Sigma(:,:,n-1)*F'+[eye(3)*randomWalk^2*dt,zeros(3);
            zeros(3),eye(3)*biasInstability^2*dt];
    end
    x(:,n) = x(:,n-1);
    
    % update
    if rem(n,5)==2
        if meaIsVec
            % vector measurement
            vPredict = zeros(3*nVecRef,1);
            H = zeros(3*nVecRef,6);
            vMeaNoise = zeros(3*nVecRef,3*nVecRef);
            for nv = 1:nVecRef
                if vecRefInertial
                    vPredict(3*(nv-1)+1:3*nv) = R(:,:,n)'*vRef(3*(nv-1)+1:3*nv);
                    H(3*(nv-1)+1:3*nv,1:3) = hat(R(:,:,n)'*vRef(3*(nv-1)+1:3*nv));
                else
                    vPredict(3*(nv-1)+1:3*nv) = R(:,:,n)*vRef(3*(nv-1)+1:3*nv);
                    H(3*(nv-1)+1:3*nv,1:3) = -R(:,:,n)*hat(vRef(3*(nv-1)+1:3*nv));
                end
                vMeaNoise(3*(nv-1)+1:3*nv,3*(nv-1)+1:3*nv) = eye(3)*meaNoise(nv);
            end
            
            K = Sigma(:,:,n)*H'*(H*Sigma(:,:,n)*H'+vMeaNoise)^-1;
            dx = K*(Mea(:,n)-vPredict);
        else
            % attitude measurement
            H = [eye(3),zeros(3)];
            if attMeaLocal
                K = Sigma(:,:,n)*H'*(H*Sigma(:,:,n)*H'+meaNoise)^-1;
            else
                K = Sigma(:,:,n)*H'*(H*Sigma(:,:,n)*H'+R(:,:,n)'*meaNoise*R(:,:,n))^-1;
            end
            dx = K*logRot(R(:,:,n)'*Mea(:,:,n),'v');
        end
        
        R(:,:,n) = R(:,:,n)*expRot(dx(1:3));
        x(:,n) = x(:,n)+dx(4:6);
        Sigma(:,:,n) = (eye(6)-K*H)*Sigma(:,:,n);
    end
    
    stepT(n-1) = toc;
end

end


%% convert distributions: Monte Carlo
function [ Sigma ] = MF2Gau( S )

N = 100000;
R = pdf_MF_sampling(S,N);

v = logRot(R,'v');
Sigma = cov(v');

end


function [ sigmaSqr ] = VM2Gau( kappa )

N = 100000;
v = pdf_VM_sampling(kappa,[0;0;1],N);

sigmaSqr = mean([var(v(1,:),1),var(v(2,:),1)]);

end


%% sampling from von Mises-Fisher distribution
function x=pdf_VM_sampling(kappa,mu,N)
% simulate von Mises-Fisher distribution on \Sph^2
% K. Mardia and P. Jupp, Directional Statistics, 2000, pp. 172

%assert(abs(norm(mu)-1)<1e-5,'keyboard');

R=[mu null(mu') ];
if det(R) < 0
    R=R(:,[1 3 2]);
end

x=zeros(3,N);
for k=1:N
   theta=inv_cdf_VMF_theta(rand,kappa);
   phi=rand*2*pi;
   
   x0=[cos(theta); sin(theta)*cos(phi); sin(theta)*sin(phi)];
   x(:,k)=R*x0;    
end

end


function theta=inv_cdf_VMF_theta(C,kappa)
theta=acos(log(exp(kappa)-2*sinh(kappa)*C)/kappa);
end

