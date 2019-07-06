function [ P,V,E,A ] = eulInt(acce, gyro, acceOld, gyroOld, Pold, Vold, Eold, dt, g)
% First order attitude/velocity/position integration using euler angles
% output: P->position, V->velocity, E->Euler angles, A->free acceleration
% input: acce->acceleration at step n;
%        gyro->angular velocity at step n;
%        acceOld->acceleration at step n-1;
%        gyroOld->angular velocity at step n-1;
%        Pold->position at step n-1;
%        Vold->velocity at step n-1;
%        Eold->Euler angle at step n-1;
%        dt->time interval
%        g->usually [0;0;-9.8].
% Reference: Eric Foxlin, Intertial Head-Tracker Sensor Fusion by a Complementary Separate-Bias Kalman Filter, 1996
% Last edited by Weixin Wang, December 18, 2018

aveAcce = 0.5*(acce+acceOld);
aveGyro = 0.5*(gyro+gyroOld);
dAcce = (acce-acceOld)/dt;
WB(1,1:3) = [1, sin(Eold(1))*tan(Eold(2)), cos(Eold(1))*tan(Eold(2))];
WB(2,1:3) = [0, cos(Eold(1)), -sin(Eold(1))];
WB(3,1:3) = [0, sin(Eold(1))/cos(Eold(2)), cos(Eold(1))/cos(Eold(2))];
VB(1,1) = cos(Eold(1))*tan(Eold(2))*gyroOld(2) - sin(Eold(1))*tan(Eold(2))*gyroOld(3);
VB(1,2) = sin(Eold(1))*cos(Eold(2))^-2*gyroOld(2) + cos(Eold(1))*cos(Eold(2))^-2*gyroOld(3);
VB(1,3) = 0;
VB(2,1) = -sin(Eold(1))*gyroOld(2)-cos(Eold(1))*gyroOld(3);
VB(2,2) = 0;
VB(2,3) = 0;
VB(3,1) = cos(Eold(1))/cos(Eold(2))*gyroOld(2) - sin(Eold(1))/cos(Eold(2))*gyroOld(3);
VB(3,2) = sin(Eold(1))*sin(Eold(2))*cos(Eold(2))^-2*gyroOld(2) + cos(Eold(1))*sin(Eold(2))*cos(Eold(2))^-2*gyroOld(3);
VB(3,3) = 0;

RM = eul2rm(Eold);
A = RM*aveAcce+g;
P = Pold + Vold*dt + 0.5*A*dt^2;
V = Vold + A*dt + 0.5*(RM*skew(aveGyro)*aveAcce+RM*dAcce)*dt^2;
E = Eold + WB*aveGyro*dt + 0.5*VB*WB*gyroOld*dt^2;


end

