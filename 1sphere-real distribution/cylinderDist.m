function [  ] = cylinderDist(  )

% parameters
miu = 0;
miu0 = pi/4;
k = 10;
rho1 = -0.5;
rho2 = 0.5;
sigma = 1;

% intermediate parameters
b1 = sigma*sqrt(k)*rho1;
b2 = sigma*sqrt(k)*rho2;
b0 = miu-b1*cos(miu0)-b2*sin(miu0);
miuc = @(theta)b0+b1*cos(theta)+b2*sin(theta);

rho = sqrt(rho1^2+rho2^2);
sigmac = sigma*sqrt(1-rho^2);

% density
f = @(theta,x)1/(2*pi*besseli(0,k))*exp(k*cos(theta-miu0))*...
    1/sqrt(2*pi*sigmac^2)*exp(-(x-miuc(theta))^2/(2*sigmac^2));

% coordinates
Nt = 500;
Nz = 500;
theta = linspace(-pi,pi,Nt);
x = cos(theta);
y = sin(theta);
z = linspace(-4,4,Nz);

% colormap
color = zeros(Nt,Nz);
for nt = 1:Nt
    for nz = 1:Nz
        color(nz,nt) = f(theta(nt),z(nz));
    end
end

% plot
figure; hold on;
surf(repmat(x,Nz,1),repmat(y,Nz,1),repmat(z',1,Nt),color,'LineStyle','none');
plot3(x,y,repmat(z(1),1,nt),'Color','k');
plot3(x,y,repmat(z(end),1,nt),'Color','k');
view([1,1,0.2]);
set(gca,'DataAspectRatio',[1,1,4]);
xlabel('cos$\theta$','Interpreter','latex');
ylabel('sin$\theta$','Interpreter','latex');

end