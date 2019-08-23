function [ miu, sigma, P, k, theta0 ] = cylinderDistOneParaMLE( x, theta, theta0 )

% default initial
if ~exist('theta0','var')
    theta0 = atan2(mean(sin(theta)),mean(cos(theta)));
    option = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-8);
    A = @(k)besseli(1,k)/besseli(0,k);
    k = fsolve(@(k)A(k)-mean(cos(theta-theta0)),1,option);
end

% step size
alpha = 0.01;

% initial log-likelihood
[miu,sigma,P,k] = fromTheta0(x,theta,theta0,k);
sigmac = sqrt(sigma^2-k*P^2);
l = -log(sigmac)-1/2/sigmac^2*mean((x-miu-k*P*sin(theta-theta0)).^2)-...
    log(besseli(0,k))+k*mean(cos(theta-theta0));

% gradient descent
i = 1;
while i==1 || abs(l-lold)>1e-10
    lold = l;
    i = i+1;
    
    % update theta0
    dtheta0 = -1/sigmac^2*mean((x-miu-k*P*sin(theta-theta0)).*...
        (k*P*cos(theta-theta0)))+k*mean(sin(theta-theta0));
    theta0 = theta0+alpha*dtheta0;
    
    % check if k<0
    Ecos = mean(cos(theta-theta0));
    if Ecos<0
        theta0 = -theta0;
    end
    
    % new log-likelihood
    [miu,sigma,P,k] = fromTheta0(x,theta,theta0,k);
    sigmac = sqrt(sigma^2-k*P^2);
    l = -log(sigmac)-1/2/sigmac^2*mean((x-miu-k*P*sin(theta-theta0)).^2)-...
        log(besseli(0,k))+k*mean(cos(theta-theta0));
    if l-lold<0
        warning('Objective function starts decreasing');
    end
end

end


function [ miu, sigma, P, k ] = fromTheta0( x, theta, theta0, k0 )

% Von Mises part
Ecos = mean(cos(theta-theta0));

option = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-8);
A = @(k)besseli(1,k)/besseli(0,k);
k = fsolve(@(k)A(k)-Ecos,k0,option);

% correlation part
Ex = mean(x);
Esin = mean(sin(theta-theta0));
covxsin = mean((x-Ex).*sin(theta-theta0));
covsinsin = mean(sin(theta-theta0).*sin(theta-theta0))-Esin^2;

P = covxsin/covsinsin/k;

% Gaussian part
covxx = mean((x-Ex).*x);

miu = Ex-k*P*Esin;
sigma = sqrt(covxx-k*P*covxsin+k*P^2);

end

