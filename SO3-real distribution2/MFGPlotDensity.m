function [  ] = MFGPlotDensity(  )

% parameters
n1 = 1;

Miu = 0;
Sigma = 1;
U = expRot([rand,0,0]);
V = expRot([0,0,0]);
S = diag([25,0,0]);
P = [0,0.7,0]/sqrt(25);

% intermediate parameters
F = U*S*V';
M = U*V';

% other intermediate parameters
Q = @(R)U'*R*V;
Miuc = @(R)Miu + P*vee(Q(R)'*S-S*Q(R));
Sigmac = Sigma-P*(trace(S)*eye(3)-S)*P';

% Normalizing constant
c1 = 1/sqrt((2*pi)^n1*det(Sigmac));
c2 = pdf_MF_normal(diag(S));

% density
f = @(x1,R)1/c1*exp(-1/2*(x1-Miuc(R))'*Sigmac^-1*(x1-Miuc(R)))*...
    1/c2*exp(trace(F*R'));

% spherical grid
Nt1 = 100;
Nt2 = 100;
theta1 = linspace(-pi,pi,Nt2);
theta2 = linspace(0,pi,Nt1);
s1 = cos(theta1)'.*sin(theta2);
s2 = sin(theta1)'.*sin(theta2);
s3 = repmat(cos(theta2),Nt1,1);

% linear grid
Nx1 = 11;
x1 = linspace(-3,3,Nx1);

% color map
c = zeros(Nt1,Nt2,Nx1,3);
for nx1 = 1:Nx1
    fRot = @(R)f(x1(nx1),R);
    c(:,:,nx1,:) = plotSO3Dist(theta1,theta2,fRot);
end

% plot Matrix Fisher
for nx1 = 1:Nx1
    figure; hold on;
    surf(s1*0.999,s2*0.999,s3*0.999,'FaceColor',[0.9,0.9,0.9],'LineStyle','none');
    colormap([linspace(1,0,100)',linspace(1,0,100)',linspace(1,1,100)';
        linspace(0,1,100)',linspace(0,1,100)',linspace(1,1,100)';
        linspace(1,1,100)',linspace(1,0,100)',linspace(1,0,100)';
        linspace(1,1,100)',linspace(0,1,100)',linspace(0,1,100)';
        linspace(1,1,100)',linspace(1,1,100)',linspace(1,0,100)']);
    alphamap(linspace(0,1,100));
    caxis([0,1]);
    alim([0,1.5]);
    
    surf(s1,s2,s3,c(:,:,nx1,1)/max(max(c(:,:,nx1,1)))*0.2,...
        'AlphaData',c(:,:,nx1,1)/max(max(c(:,:,nx1,1)))*1.5,...
        'FaceAlpha','flat','AlphaDataMapping','scaled','LineStyle','none');
    surf(s1,s2,s3,c(:,:,nx1,2)/max(max(c(:,:,nx1,2)))*0.2+0.4,...
        'AlphaData',c(:,:,nx1,2)/max(max(c(:,:,nx1,2)))*1.25,...
        'FaceAlpha','flat','AlphaDataMapping','scaled','LineStyle','none');
    surf(s1,s2,s3,c(:,:,nx1,3)/max(max(c(:,:,nx1,3)))*0.2+0.8,...
        'AlphaData',c(:,:,nx1,3)/max(max(c(:,:,nx1,3))),...
        'FaceAlpha','flat','AlphaDataMapping','scaled','LineStyle','none');
    
    plot3([0,M(1,1)*1.1],[0,M(2,1)*1.1],[0,M(3,1)*1.1],'b');
    plot3([0,M(1,2)*1.1],[0,M(2,2)*1.1],[0,M(3,2)*1.1],'r');
    plot3([0,M(1,3)*1.1],[0,M(2,3)*1.1],[0,M(3,3)*1.1],'y');
    
    axis equal;
    view([1,1,1]);
end

% % linear grid
% Nx1 = 1001;
% x1 = linspace(-4,4,Nx1);
% 
% % rotation grid
% Nt1 = 11;
% theta1 = linspace(-pi/8,pi/8,Nt1);
% 
% % plot Gaussian
% M = U*V';
% 
% fLinear = zeros(Nx1,Nt1,3);
% MiuLinear = zeros(Nt1,3);
% for i = 1:3
%     n = zeros(3,1);
%     for nt1 = 1:Nt1
%         n(i) = theta1(nt1);
%         MiuLinear(nt1,i) = Miuc(M*expRot(V*n));
%         for nx1 = 1:Nx1
%             fLinear(nx1,nt1,i) = f(x1(nx1),M*expRot(V*n));
%         end
%     end
% end
% fmax = max(max(max(fLinear)));
% 
% for i = 1:3
%     for nt1 = 1:Nt1
%         figure; hold on;
%         plot(x1,fLinear(:,nt1,i));
%         
%         xlim([-4,4]);
%         ylim([0,fmax]);
%         
%         str = strcat('$\mu=$',num2str(MiuLinear(nt1,i)),', $\sigma=$',num2str(Sigmac));
%         annotation('textbox','String',str,'Interpreter','latex');
%     end
% end

end

