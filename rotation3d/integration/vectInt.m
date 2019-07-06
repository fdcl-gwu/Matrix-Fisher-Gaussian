function [ P,V,RV,A ] = vectInt( acce, gyro, acceOld, gyroOld, Pold, Vold, RVold, dt, g )
% First order attitude/velocity/position integration using rotation vector
% still under development

RM = expRM(RVold);
aveAcce = 0.5*(acce+acceOld);
aveGyro = 0.5*(gyro+gyroOld);
dAcce = (acce-acceOld)/dt;
A = RM*aveAcce+g;
P = Pold + Vold*dt + 0.5*A*dt^2;
V = Vold + A*dt + 0.5*(RM*skew(aveGyro)*aveAcce+RM*dAcce)*dt^2;
RV = RVold + (aveGyro + 0.5*cross(RVold,aveGyro) + 1/norm(RVold)^2*(1-RVold/2*1/tan(RVold/2))*cross(RVold,cross(RVold,aveGyro)))*dt;

end

