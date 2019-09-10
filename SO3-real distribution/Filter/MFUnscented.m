function [ R, F ] = MFUnscented( gyro, RInit, RMea, sf )

filePath = mfilename('fullpath');
pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    addpath(getAbsPath('..\..\rotation3d',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end

dt = 1/sf;
N = size(gyro,2);

% noise parameters
randomWalk = 10*pi/180;
rotMeaNoise = 0.2;

SM = Gau2MF(rotMeaNoise);

% data containers
R = zeros(3,3,N);
s = zeros(3,N);

% initialize
R(:,:,1) = RInit*expRot([pi,0,0]);
s(:,1) = [0;0;0];
F = RInit*diag(s(:,1));

% filter iterations
for n = 2:N
    % unscented transform from last step
    [Rl,wl] = MFGetSigmaPoints(F);
    [xav,wav] = GGetSigmaPoints([0;0;0],eye(3)*randomWalk^2*sf);
    
    % propagate sigma points
    Rp = zeros(3,3,7*7);
    wp = zeros(1,7*7);
    for i = 1:7
        for j = 1:7
            ind = 7*(i-1)+j;
            av = (gyro(:,n-1)+gyro(:,n))/2;
            Rp(:,:,ind) = Rl(:,:,i)*expRot((av-xav(:,j))*dt);
            wp(ind) = wl(i)*wav(j);
        end
    end
    
    % recover prior distribution
    ER = sum(Rp.*reshape(wp,1,1,[]),3);
    [U,D,V] = psvd(ER);
    S = diag(pdf_MF_M2S(diag(D)));
    F = U*S*V';
    
    % update
    if rem(n,5)==0
        F = F+RMea(:,:,n)*SM;
    end
    
    [U,S,V] = psvd(F);
    R(:,:,n) = U*V';
    s(:,n) = diag(S);
end

if ~any(strcmp(pathCell,getAbsPath('..\..\rotation3d',filePath)))
    rmpath(getAbsPath('..\..\rotation3d',filePath));
end
if ~any(strcmp(pathCell,getAbsPath('..\Matrix-Fisher-Distribution',filePath)))
    addpath(getAbsPath('..\Matrix-Fisher-Distribution',filePath));
end

end


function [ S ] = Gau2MF( sigma )

N = 100000;
v = sigma*randn(3,N);

R = expRot(v);
ER = mean(R,3);
[~,D,~] = psvd(ER);

S = pdf_MF_M2S(diag(D));
S = eye(3)*mean(S);

end

