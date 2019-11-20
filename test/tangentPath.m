addpath('..\rotation3d');

% grid
N1 = 100;
p1 = [linspace(0,pi/2,N1)',zeros(N1,2)];
N2 = 400;
p2 = [ones(N2,1)*pi/2,linspace(0,pi*2,N2)',zeros(N2,1)];

RInit = eye(3);
R = zeros(3,3,N1+N2);
for n1 = 1:N1
    R(:,:,n1) = RInit*expRot(p1(n1,:));
end
for n2 = 1:N2
    R(:,:,N1+n2) = RInit*expRot(p2(n2,:));
end

figure; hold on;
view([1,1,1]);
for n1 = 1:N1
    plot3([0,R(1,1,n1)],[0,R(2,1,n1)],[0,R(3,1,n1)],'b');
    plot3([0,R(1,2,n1)],[0,R(2,2,n1)],[0,R(3,2,n1)],'r');
    plot3([0,R(1,3,n1)],[0,R(2,3,n1)],[0,R(3,3,n1)],'y');
    pause(0.1);
    axis equal;
    drawnow;
end
for n2 = 1:N2
    plot3([0,R(1,1,N1+n2)],[0,R(2,1,N1+n2)],[0,R(3,1,N1+n2)],'b');
    plot3([0,R(1,2,N1+n2)],[0,R(2,2,N1+n2)],[0,R(3,2,N1+n2)],'r');
    plot3([0,R(1,3,N1+n2)],[0,R(2,3,N1+n2)],[0,R(3,3,N1+n2)],'y');
    pause(0.1);
    axis equal;
    drawnow;
end

rmpath('..\rotation3d');
