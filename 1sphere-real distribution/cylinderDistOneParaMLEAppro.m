function [ miu, sigma, P, k, theta0 ] = cylinderDistOneParaMLEAppro( x, theta )

% Von Mises part
theta0 = atan2(mean(sin(theta)),mean(cos(theta)));
option = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-8);
A = @(k)besseli(1,k)/besseli(0,k);
k = fsolve(@(k)A(k)-mean(cos(theta-theta0)),1,option);

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

