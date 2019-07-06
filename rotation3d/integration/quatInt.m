function [ P,V,Q,A ] = quatInt(acce, gyro, acceOld, gyroOld, Pold, Vold, Qold, dt, g)
% First order attitude/velocity/position integration using quaternions
% output: P->position, V->velocity, Q->quaternion, A->free acceleration
% input: acce->acceleration at step n;
%        gyro->angular velocity at step n;
%        acceOld->acceleration at step n-1;
%        gyroOld->angular velocity at step n-1;
%        Pold->position at step n-1;
%        Vold->velocity at step n-1;
%        Qold->quaternion at step n-1;
%        dt->time interval
%        g->usually [0;0;-9.8].
% Reference: Joan Sola, Quaternion kinematics for the error-state Kalman filter, 2017
% Last edited by Weixin Wang, December 18, 2018

RM = QuatToRM(Qold);
aveAcce = 0.5*(acce+acceOld);
aveGyro = 0.5*(gyro+gyroOld);
A = RM*aveAcce+g;
P = Pold + Vold*dt + 0.5*A*dt^2;
V = Vold + A*dt + 0.5*RM*skew(aveGyro)*(RM'*A)*dt^2;
Q = QuatMultiply(Qold,expQuat(aveGyro*dt)+dt^2/24*[0;cross(gyroOld,gyro)]);

end

