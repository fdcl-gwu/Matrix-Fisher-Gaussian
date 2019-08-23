function [ ] = cylinderMLETest( x, theta )

Nt0 = 10000;
theta0 = linspace(-pi,pi,Nt0);

l = zeros(Nt0,1);
kAll = zeros(Nt0,1);
parfor i = 1:Nt0
    [miu,sigma,P,k] = fromTheta0(x,theta,theta0(i));
    sigmac = sqrt(sigma^2-k*P^2);
    l(i) = -log(sigmac)-1/2/sigmac^2*mean((x-miu-k*P*sin(theta-theta0(i))).^2)-...
        log(besseli(0,k))+k*mean(cos(theta-theta0(i)));
    kAll(i) = k;
end

theta0VM = atan2(mean(sin(theta)),mean(cos(theta)));

figure; hold on;
plot(theta0,l);
plot(theta0,kAll);
yLim = get(gca,'YLim');
plot([theta0VM,theta0VM],yLim);

zeroInd = find(kAll(1:end-1).*kAll(2:end)<0);
for i = 1:length(zeroInd)
    plot([theta0(zeroInd(i)),theta0(zeroInd(i))],yLim,'k');
end

end


function [ miu, sigma, P, k ] = fromTheta0( x, theta, theta0 )

% Von Mises part
Ecos = mean(cos(theta-theta0));

option = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-8);
A = @(k)besseli(1,k)/besseli(0,k);
k = fsolve(@(k)A(k)-Ecos,1,option);

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

