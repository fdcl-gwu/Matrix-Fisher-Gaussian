function [ x, theta ] = cylinderSampling( miu, sigma, P, k, theta0, Ns )

% Von Mises part
theta = VMSampling(k,theta0,Ns);

% Gaussian part
if size(P,2) == 1
    sigmac = sqrt(sigma^2-k*P^2);
    miuc = miu+k*P*sin(theta-theta0);
elseif size(P,2) == 2
    sigmac = sqrt(sigma^2-k*(P*P'));
    miuc = miu+k*P*[cos(theta)-cos(theta0);sin(theta)-sin(theta0)];
end

x = randn(1,Ns).*sigmac+miuc;

end

