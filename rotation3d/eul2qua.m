% Author: Lauro Ojeda, 2003-2015
% Bibliography:
% Collinson: Introduction to Avionics Systems, 2nd Edition
% Collinson pp 266
function Q = eul2qua(Att)

if size(Att,2) ~= 3
    Att = Att';
end

phi=Att(:,1);
the=Att(:,2);
psi=Att(:,3);
cos_phi=cos(phi/2.0);
sin_phi=sin(phi/2.0);
cos_the=cos(the/2.0);
sin_the=sin(the/2.0);
cos_psi=cos(psi/2.0);
sin_psi=sin(psi/2.0);
a=cos_phi.*cos_the.*cos_psi+sin_phi.*sin_the.*sin_psi;
b=sin_phi.*cos_the.*cos_psi-cos_phi.*sin_the.*sin_psi;
c=cos_phi.*sin_the.*cos_psi+sin_phi.*cos_the.*sin_psi;
d=cos_phi.*cos_the.*sin_psi-sin_phi.*sin_the.*cos_psi;
Q=[a,b,c,d];

end
