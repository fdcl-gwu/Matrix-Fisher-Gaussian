function [ Miu, Sigma, P, U, S, V ] = SO3RealMLE( x, R, w, U, S, V )
% let x be N-by-Ns, R be 3-by-3-by-Ns

pathCell = regexp(path, pathsep, 'split');
if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    addpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    addpath('..\rotation3d');
end

N = size(x,1);
Ns = size(R,3);

if ~exist('w','var') || isempty(w)
    w = ones(1,Ns)/Ns;
end

% default initial estimates
if ~exist('U0','var')
%     ER = sum(R.*reshape(w,1,1,Ns),3);
%     [U,D,V] = usvd(ER);
%     S = diag(pdf_MF_M2S(diag(D)));
    U = eye(3);
    V = eye(3);
    S = eye(3);
end

% step size
kS = 1;
kUV = 0.001;

% initial log likelihood
[Miu,Sigma,P,fR,Sigmac,Q] = fromF(x,R,w,U,S,V);
ER = sum(R.*reshape(w,1,1,Ns),3);
l = -log(pdf_MF_normal(diag(S))) + trace(U*S*V'*ER.') -...
    1/2*log(det(Sigmac));
parfor ns = 1:Ns
    l = l-1/2*w(ns)*(x(:,ns)-Miu-P*fR(:,ns))'*Sigmac^-1*...
        (x(:,ns)-Miu-P*fR(:,ns));
end
l

% gradient descent iteration
i = 1;
while i==1 || abs(l-loldS)>1e-6
    i = i+1;
    loldS = l;
    
    j = 1;
    while j==1 || abs(l-loldU)>1e-6
        j = j+1;
        loldU = l;
        
        k = 1;
        while k==1 || abs(l-loldV)>1e-8
            k = k+1;
            loldV = l;
            
            % gradient descent
            dS = [sum(Q(1,1,:).*reshape(w,1,1,Ns));
                sum(Q(2,2,:).*reshape(w,1,1,Ns));
                sum(Q(3,3,:).*reshape(w,1,1,Ns))]...
                -pdf_MF_normal_deriv(diag(S))/pdf_MF_normal(diag(S));
            parfor ns = 1:Ns
                dS = dS+w(ns)*[0, -Q(3,1,ns), Q(2,1,ns);
                    Q(3,2,ns), 0, -Q(1,2,ns);
                    -Q(2,3,ns), Q(1,3,ns), 0]*...
                    P'*Sigmac^-1*(x(:,ns)-Miu-P*fR(:,ns));
            end
            S = S+kS*diag(dS);

            % new log-likelihood
            [Miu,Sigma,P,fR,Sigmac,Q] = fromF(x,R,w,U,S,V);
            ER = sum(R.*reshape(w,1,1,Ns),3);
            l = -log(pdf_MF_normal(diag(S))) + trace(U*S*V'*ER.') -...
                1/2*log(det(Sigmac));
            parfor ns = 1:Ns
                l = l-1/2*w(ns)*(x(:,ns)-Miu-P*fR(:,ns))'*Sigmac^-1*...
                    (x(:,ns)-Miu-P*fR(:,ns));
            end
            l
        end
        
        % gradient descent
        EQ = sum(Q.*reshape(w,1,1,Ns),3);
        dU = vee(EQ*S-S*EQ');
        parfor ns = 1:Ns
            dU = dU+w(ns)*[-Q(2,2,ns)*S(2,2)-Q(3,3,ns)*S(3,3), Q(2,1,ns)*S(1,1), Q(3,1,ns)*S(1,1);
                Q(1,2,ns)*S(2,2), -Q(1,1,ns)*S(1,1)-Q(3,3,ns)*S(1,1), Q(3,2,ns)*S(2,2);
                Q(1,3,ns)*S(3,3), Q(2,3,ns)*S(3,3), -Q(1,1,ns)*S(1,1)-Q(2,2,ns)*S(2,2)]*...
                P'*Sigmac^-1*(x(:,ns)-Miu-P*fR(:,ns));
        end
        U = U*expRM(kUV*dU);

        % new log-likelihood
        [Miu,Sigma,P,fR,Sigmac,Q] = fromF(x,R,w,U,S,V);
        ER = sum(R.*reshape(w,1,1,Ns),3);
        l = -log(pdf_MF_normal(diag(S))) + trace(U*S*V'*ER.') -...
            1/2*log(det(Sigmac));
        parfor ns = 1:Ns
            l = l-1/2*w(ns)*(x(:,ns)-Miu-P*fR(:,ns))'*Sigmac^-1*...
                (x(:,ns)-Miu-P*fR(:,ns));
        end
        l
    end
    
    % gradient descent
    EQ = sum(Q.*reshape(w,1,1,Ns),3);
    dV = vee(EQ'*S-S*EQ);
    parfor ns = 1:Ns
        dV = dV+w(ns)*[Q(2,2,ns)*S(3,3)+Q(3,3,ns)*S(2,2), -Q(1,2,ns)*S(3,3), -Q(1,3,ns)*S(2,2);
            -Q(2,1,ns)*S(3,3), Q(1,1,ns)*S(3,3)+Q(3,3,ns)*S(1,1), -Q(2,3,ns)*S(1,1);
            -Q(3,1,ns)*S(2,2), -Q(3,2,ns)*S(1,1), Q(1,1,ns)*S(2,2)+Q(2,2,ns)*S(1,1)]*...
            P'*Sigmac^-1*(x(:,ns)-Miu-P*fR(:,ns));
    end
    V = V*expRM(kUV*dV);

    % new log-likelihood
    [Miu,Sigma,P,fR,Sigmac,Q] = fromF(x,R,w,U,S,V);
    ER = sum(R.*reshape(w,1,1,Ns),3);
    l = -log(pdf_MF_normal(diag(S))) + trace(U*S*V'*ER.') -...
        1/2*log(det(Sigmac));
    parfor ns = 1:Ns
        l = l-1/2*w(ns)*(x(:,ns)-Miu-P*fR(:,ns))'*Sigmac^-1*...
            (x(:,ns)-Miu-P*fR(:,ns));
    end
    l
end

if ~any(strcmp(pathCell,getAbsPath('Matrix-Fisher-Distribution')))
    rmpath('Matrix-Fisher-Distribution');
end
if ~any(strcmp(pathCell,getAbsPath('..\rotation3d')))
    rmpath('..\rotation3d');
end

end


function [ Miu, Sigma, P, fR, Sigmac, Q ] = fromF( x, R, w, U, S, V )

N = size(x,1);
Ns = size(R,3);

fR = zeros(3,Ns);
Q = zeros(3,3,Ns);
parfor ns = 1:Ns
    Q(:,:,ns) = U.'*R(:,:,ns)*V;
    fR(:,ns) = vee(Q(:,:,ns)*S-S*Q(:,:,ns).')/sqrt(2);
end

% correlation part
Ex = sum(x.*w,2);
EfR = sum(fR.*w,2);
covxfR = sum(reshape(x,N,1,Ns).*reshape(fR,1,3,Ns).*...
    reshape(w,1,1,Ns),3)-Ex*EfR.';
covfRfR = sum(reshape(fR,3,1,Ns).*reshape(fR,1,3,Ns).*...
    reshape(w,1,1,Ns),3)-EfR*EfR.';
P = covxfR*covfRfR^-1;

% Gaussian part
covxx = sum(reshape(x,N,1,Ns).*reshape(x,1,N,Ns).*reshape(w,1,1,Ns),3)-...
    Ex*Ex.';
Sigma2Inv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)])/2;
Miu = Ex-P*EfR;
Sigmac = covxx-P*covxfR';
Sigma = Sigmac+P*Sigma2Inv*P';

end

