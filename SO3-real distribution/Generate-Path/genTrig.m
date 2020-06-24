function [ gyroMea, Mea, RTrue, biasTrue ] = genTrig( t, sf, parameters )

time = (0:1/sf:t);
N = length(time);

%% parameters
% motion parameters
E.fr = 0.35; E.fp = 0.35; E.fh = 0.35;
E.magr = pi; E.magp = pi/2; E.magh = pi;

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

%% true attitude
roll = @(t)E.magr*sin(E.fr*2*pi*t);
pitch = @(t)E.magp*sin(E.fp*2*pi*t);
yaw = @(t)E.magh*sin(E.fh*2*pi*t);

sr2 = @(t)sin(roll(t)/2); sp2 = @(t)sin(pitch(t)/2); sy2 = @(t)sin(yaw(t)/2);
cr2 = @(t)cos(roll(t)/2); cp2 = @(t)cos(pitch(t)/2); cy2 = @(t)cos(yaw(t)/2);
q1 = @(t)cr2(t).*cp2(t).*cy2(t)+sr2(t).*sp2(t).*sy2(t);
q2 = @(t)sr2(t).*cp2(t).*cy2(t)-cr2(t).*sp2(t).*sy2(t);
q3 = @(t)cr2(t).*sp2(t).*cy2(t)+sr2(t).*cp2(t).*sy2(t);
q4 = @(t)cr2(t).*cp2(t).*sy2(t)-sr2(t).*sp2(t).*cy2(t);
qTrue = [q1(time);q2(time);q3(time);q4(time)];
RTrue = qua2rot(qTrue);

%% true angular velocity
dr = @(t)E.magr*E.fr*2*pi*cos(E.fr*2*pi*t);
dp = @(t)E.magp*E.fp*2*pi*cos(E.fp*2*pi*t);
dy = @(t)E.magh*E.fh*2*pi*cos(E.fh*2*pi*t);

if omegaLocal
    omega = @(t) [dr(t)-sin(pitch(t)).*dy(t);
        cos(roll(t)).*dp(t)+sin(roll(t)).*cos(pitch(t)).*dy(t);
        -sin(roll(t)).*dp(t)+cos(roll(t)).*cos(pitch(t)).*dy(t)];
    gyro = omega(time);
else
    omega = @(t) [cos(yaw(t)).*cos(pitch(t)).*dr(t) - sin(yaw(t)).*dp(t);
        sin(yaw(t)).*cos(pitch(t)).*dr(t) + cos(yaw(t)).*dp(t);
        -sin(pitch(t)).*dr(t) + dy(t)];
    gyro = omega(time);
end

%% add noise
biasNoise = randn(3,N)*biasInstability*sqrt(sf);
biasTrue = cumsum(biasNoise/sf,2);

gyroNoise = randn(3,N)*randomWalk*sqrt(sf);
gyroMea = gyro+gyroNoise+biasTrue;

%% measurement
if meaIsVec
    Mea = zeros(3*nVecRef,N);
    if GaussMea
        vecNoise = zeros(3*nVecRef,N);
        for nv = 1:nVecRef
            vecNoise(3*(nv-1)+1:3*nv,:) = randn(3,N)*sqrt(meaNoise(nv));
        end
        
        for n = 1:N
            for nv = 1:nVecRef
                if vecRefInertial
                    Mea(3*(nv-1)+1:3*nv,n) = RTrue(:,:,n)'*vRef(3*(nv-1)+1:3*nv)...
                        + vecNoise(3*(nv-1)+1:3*nv,n);
                else
                    Mea(3*(nv-1)+1:3*nv,n) = RTrue(:,:,n)*vRef(3*(nv-1)+1:3*nv)...
                        + vecNoise(3*(nv-1)+1:3*nv,n);
                end
                Mea(3*(nv-1)+1:3*nv,n) = Mea(3*(nv-1)+1:3*nv,n)./sqrt(sum(Mea(3*(nv-1)+1:3*nv,n).^2));
            end
        end
    else
        for n = 1:N
            for nv = 1:nVecRef
                if vecRefInertial
                    Mea(3*(nv-1)+1:3*nv,n) = pdf_VM_sampling(meaNoise(nv),RTrue(:,:,n)'*vRef(3*(nv-1)+1:3*nv),1);
                else
                    Mea(3*(nv-1)+1:3*nv,n) = pdf_VM_sampling(meaNoise(nv),RTrue(:,:,n)*vRef(3*(nv-1)+1:3*nv),1);
                end
            end
        end
    end
else
    if GaussMea
        RNoise = expRot(mvnrnd([0;0;0],meaNoise,N));
    else
        RNoise = pdf_MF_sampling(meaNoise,N);
    end
    
    if attMeaLocal
        Mea = mulRot(RTrue,RNoise);
    else
        Mea = mulRot(RNoise,RTrue);
    end
end

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

