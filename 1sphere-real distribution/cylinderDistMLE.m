function [ miu, sigma, P, k, theta0 ] = cylinderDistMLE( x, theta )

% Von Mises part
Ecos = mean(cos(theta));
Esin = mean(sin(theta));

theta0 = atan2(Esin,Ecos);

option = optimoptions('fsolve','Display','off');
A = @(k)besseli(1,k)/besseli(0,k);
k = fsolve(@(k)A(k)-sqrt(Ecos^2+Esin^2),1,option);

% correlation part
Ex = mean(x);
covxcos = mean((x-Ex).*(cos(theta)-Ecos));
covxsin = mean((x-Ex).*(sin(theta)-Esin));
covcoscos = mean((cos(theta)-Ecos).*(cos(theta)-Ecos));
covsinsin = mean((sin(theta)-Esin).*(sin(theta)-Esin));
covcossin = mean((cos(theta)-Ecos).*(sin(theta)-Esin));

P = [covxcos,covxsin]*[covcoscos,covcossin;covcossin,covsinsin]^-1/k;

% Gaussian part
covxx = mean((x-Ex).*(x-Ex));

miu = Ex-P*[Ecos;Esin];
sigma = sqrt(covxx-P*k*[covxcos;covxsin]+k*(P*P'));

end

