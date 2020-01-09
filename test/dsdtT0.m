function [ output ] = dsdtT0( inds, indT )

% initialize symbolic variables
syms U11 U12 U13 U21 U22 U23 U31 U32 U33;
U = [U11,U12,U13;U21,U22,U23;U31,U32,U33];

syms V11 V12 V13 V21 V22 V23 V31 V32 V33;
V = [V11,V12,V13;V21,V22,V23;V31,V32,V33];

syms s1 s2 s3;
s = [s1,s2,s3];

OmegaU = sym(zeros(3,3,3,3));
for i = 1:3
    for j = 1:3
        OmegaU(i,j,1,2) = sym(strcat('OmegaU',num2str(i),num2str(j),'12'));
        OmegaU(i,j,1,3) = sym(strcat('OmegaU',num2str(i),num2str(j),'13'));
        OmegaU(i,j,2,3) = sym(strcat('OmegaU',num2str(i),num2str(j),'23'));
        OmegaU(i,j,2,1) = -OmegaU(i,j,1,2);
        OmegaU(i,j,3,1) = -OmegaU(i,j,1,3);
        OmegaU(i,j,3,2) = -OmegaU(i,j,2,3);
    end
end

OmegaV = sym(zeros(3,3,3,3));
for i = 1:3
    for j = 1:3
        OmegaV(i,j,1,2) = sym(strcat('OmegaV',num2str(i),num2str(j),'12'));
        OmegaV(i,j,1,3) = sym(strcat('OmegaV',num2str(i),num2str(j),'13'));
        OmegaV(i,j,2,3) = sym(strcat('OmegaV',num2str(i),num2str(j),'23'));
        OmegaV(i,j,2,1) = -OmegaV(i,j,1,2);
        OmegaV(i,j,3,1) = -OmegaV(i,j,1,3);
        OmegaV(i,j,3,2) = -OmegaV(i,j,2,3);
    end
end

% calculate
output = dsdt(inds,indT,U,V,s,OmegaU,OmegaV);

% solve for OmegaU and OmegaV
for i = 1:3
    for j = 1:3
        temp12 = linsolve([s2,s1;s1,s2],[U(i,1)*V(j,2);-U(i,2)*V(j,1)]);
        temp13 = linsolve([s3,s1;s1,s3],[U(i,1)*V(j,3);-U(i,3)*V(j,1)]);
        temp23 = linsolve([s3,s2;s2,s3],[U(i,2)*V(j,3);-U(i,3)*V(j,2)]);
        
        OmegaU(i,j,1,2) = temp12(1);
        OmegaU(i,j,1,3) = temp13(1);
        OmegaU(i,j,2,3) = temp23(1);
        OmegaU(i,j,2,1) = -temp12(1);
        OmegaU(i,j,3,1) = -temp13(1);
        OmegaU(i,j,3,2) = -temp23(1);
        
        OmegaV(i,j,1,2) = temp12(2);
        OmegaV(i,j,1,3) = temp13(2);
        OmegaV(i,j,2,3) = temp23(2);
        OmegaV(i,j,2,1) = -temp12(2);
        OmegaV(i,j,3,1) = -temp13(2);
        OmegaV(i,j,3,2) = -temp23(2);
    end
end

% evaluate at T=0
for i = 1:3
    for j = 1:3
        output = subs(output,str2sym(strcat('OmegaU',num2str(i),num2str(j),'12')),OmegaU(i,j,1,2));
        output = subs(output,str2sym(strcat('OmegaU',num2str(i),num2str(j),'13')),OmegaU(i,j,1,3));
        output = subs(output,str2sym(strcat('OmegaU',num2str(i),num2str(j),'23')),OmegaU(i,j,2,3));
        output = subs(output,str2sym(strcat('OmegaU',num2str(i),num2str(j),'21')),OmegaU(i,j,2,1));
        output = subs(output,str2sym(strcat('OmegaU',num2str(i),num2str(j),'31')),OmegaU(i,j,3,1));
        output = subs(output,str2sym(strcat('OmegaU',num2str(i),num2str(j),'32')),OmegaU(i,j,3,2));
        
        output = subs(output,str2sym(strcat('OmegaV',num2str(i),num2str(j),'12')),OmegaV(i,j,1,2));
        output = subs(output,str2sym(strcat('OmegaV',num2str(i),num2str(j),'13')),OmegaV(i,j,1,3));
        output = subs(output,str2sym(strcat('OmegaV',num2str(i),num2str(j),'23')),OmegaV(i,j,2,3));
        output = subs(output,str2sym(strcat('OmegaV',num2str(i),num2str(j),'21')),OmegaV(i,j,2,1));
        output = subs(output,str2sym(strcat('OmegaV',num2str(i),num2str(j),'31')),OmegaV(i,j,3,1));
        output = subs(output,str2sym(strcat('OmegaV',num2str(i),num2str(j),'32')),OmegaV(i,j,3,2));
    end
end

output = subs(output,[U11,U12,U13;U21,U22,U23;U31,U32,U33],eye(3));
output = subs(output,[V11,V12,V13;V21,V22,V23;V31,V32,V33],eye(3));

end


function [ output ] = dsdt( inds, indT, U, V, s, OmegaU, OmegaV )

if size(indT,1)==1
    output = U(indT(1,1),inds)*V(indT(1,2),inds);
else
    N = size(indT(2:end,:),1);
    output = U(indT(1,1),inds)*dVdt([indT(1,2),inds],indT(2:end,:),U,V,s,OmegaU,OmegaV);
    
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dUInd = comb(nc,:)+1;
            dVInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                output = output + dUdt([indT(1,1),inds],indT(2:end,:),U,V,s,OmegaU,OmegaV)...
                    *V(indT(1,2),inds);
            else
                output = output + dUdt([indT(1,1),inds],indT(dUInd,:),U,V,s,OmegaU,OmegaV)...
                    *dVdt([indT(1,2),inds],indT(dVInd,:),U,V,s,OmegaU,OmegaV);
            end
        end
    end
end

end


function [ output ] = dUdt( indU, indT, U, V, s, OmegaU, OmegaV )

ext = setdiff([1,2,3],indU(2));

if size(indT,1)==1
    output = U(indU(1),ext(1))*OmegaU(indT(1,1),indT(1,2),ext(1),indU(2))...
        + U(indU(1),ext(2))*OmegaU(indT(1,1),indT(1,2),ext(2),indU(2));
else
    N = size(indT(2:end,:),1);
    
    % DUOmegaU1
    DUOmegaU1 = U(indU(1),ext(1))*dOmegaUdt(indT(1,:),[ext(1),indU(2)],indT(2:end,:),U,V,s,OmegaU,OmegaV);
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dUInd = comb(nc,:)+1;
            dOmegaUInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DUOmegaU1 = DUOmegaU1 + dUdt([indU(1),ext(1)],indT(2:end,:),U,V,s,OmegaU,OmegaV)...
                    *OmegaU(indT(1,1),indT(1,2),ext(1),indU(2));
            else
                DUOmegaU1 = DUOmegaU1 + dUdt([indU(1),ext(1)],indT(dUInd,:),U,V,s,OmegaU,OmegaV)...
                    *dOmegaUdt(indT(1,:),[ext(1),indU(2)],indT(dOmegaUInd,:),U,V,s,OmegaU,OmegaV);
            end
        end
    end
    
    % DUOmegaU2
    DUOmegaU2 = U(indU(1),ext(2))*dOmegaUdt(indT(1,:),[ext(2),indU(2)],indT(2:end,:),U,V,s,OmegaU,OmegaV);
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dUInd = comb(nc,:)+1;
            dOmegaUInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DUOmegaU2 = DUOmegaU2 + dUdt([indU(1),ext(2)],indT(2:end,:),U,V,s,OmegaU,OmegaV)...
                    *OmegaU(indT(1,1),indT(1,2),ext(2),indU(2));
            else
                DUOmegaU2 = DUOmegaU2 + dUdt([indU(1),ext(2)],indT(dUInd,:),U,V,s,OmegaU,OmegaV)...
                    *dOmegaUdt(indT(1,:),[ext(2),indU(2)],indT(dOmegaUInd,:),U,V,s,OmegaU,OmegaV);
            end
        end
    end
    
    output = DUOmegaU1+DUOmegaU2;
end

end


function [ output ] = dVdt( indV, indT, U, V, s, OmegaU, OmegaV )

ext = setdiff([1,2,3],indV(2));

if size(indT,1)==1
    output = -V(indV(1),ext(1))*OmegaV(indT(1,1),indT(1,2),ext(1),indV(2))...
        - V(indV(1),ext(2))*OmegaV(indT(1,1),indT(1,2),ext(2),indV(2));
else
    N = size(indT(2:end,:),1);
    
    % DVOmegaU1
    DVOmegaV1 = V(indV(1),ext(1))*dOmegaVdt(indT(1,:),[ext(1),indV(2)],indT(2:end,:),U,V,s,OmegaU,OmegaV);
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dVInd = comb(nc,:)+1;
            dOmegaVInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DVOmegaV1 = DVOmegaV1 + dVdt([indV(1),ext(1)],indT(2:end,:),U,V,s,OmegaU,OmegaV)...
                    *OmegaV(indT(1,1),indT(1,2),ext(1),indV(2));
            else
                DVOmegaV1 = DVOmegaV1 + dVdt([indV(1),ext(1)],indT(dVInd,:),U,V,s,OmegaU,OmegaV)...
                    *dOmegaVdt(indT(1,:),[ext(1),indV(2)],indT(dOmegaVInd,:),U,V,s,OmegaU,OmegaV);
            end
        end
    end
    
    % DVOmegaV2
    DVOmegaV2 = V(indV(1),ext(2))*dOmegaVdt(indT(1,:),[ext(2),indV(2)],indT(2:end,:),U,V,s,OmegaU,OmegaV);
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dVInd = comb(nc,:)+1;
            dOmegaVInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DVOmegaV2 = DVOmegaV2 + dVdt([indV(1),ext(2)],indT(2:end,:),U,V,s,OmegaU,OmegaV)...
                    *OmegaV(indT(1,1),indT(1,2),ext(2),indV(2));
            else
                DVOmegaV2 = DVOmegaV2 + dVdt([indV(1),ext(2)],indT(dVInd,:),U,V,s,OmegaU,OmegaV)...
                    *dOmegaVdt(indT(1,:),[ext(2),indV(2)],indT(dOmegaVInd,:),U,V,s,OmegaU,OmegaV);
            end
        end
    end
    
    output = -DVOmegaV1-DVOmegaV2;
end

end


function [ output ] = dOmegaUdt( indOsup, indOsub, indT, U, V, s, OmegaU, OmegaV )

N = size(indT,1);

% DslOmegaUkl
DslOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaUkl = DslOmegaUkl + dsdt(indOsub(2),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaU(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DslOmegaUkl = DslOmegaUkl + dsdt(indOsub(2),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaUdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DskOmegaVkl
DskOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaVkl = DskOmegaVkl + dsdt(indOsub(1),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaV(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DskOmegaVkl = DskOmegaVkl + dsdt(indOsub(1),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaVdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

%DskOmegaUkl
DskOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaUkl = DskOmegaUkl + dsdt(indOsub(1),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaU(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DskOmegaUkl = DskOmegaUkl + dsdt(indOsub(1),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaUdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DslOmegaVkl
DslOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaVkl = DslOmegaVkl + dsdt(indOsub(2),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaV(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DslOmegaVkl = DslOmegaVkl + dsdt(indOsub(2),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaVdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DUikVjl
DUikVjl = U(indOsup(1),indOsub(1))*...
            dVdt([indOsup(2),indOsub(2)],indT,U,V,s,OmegaU,OmegaV);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUikVjl = DUikVjl + dUdt([indOsup(1),indOsub(1)],indT,U,V,s,OmegaU,OmegaV)*...
                V(indOsup(2),indOsub(2));
        else
            DUikVjl = DUikVjl + dUdt([indOsup(1),indOsub(1)],indT(dUInd,:),U,V,s,OmegaU,OmegaV)*...
                dVdt([indOsup(2),indOsub(2)],indT(dVInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DUilVjk
DUilVjk = U(indOsup(1),indOsub(2))*...
            dVdt([indOsup(2),indOsub(1)],indT,U,V,s,OmegaU,OmegaV);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUilVjk = DUilVjk + dUdt([indOsup(1),indOsub(2)],indT,U,V,s,OmegaU,OmegaV)*...
                V(indOsup(2),indOsub(1));
        else
            DUilVjk = DUilVjk + dUdt([indOsup(1),indOsub(2)],indT(dUInd,:),U,V,s,OmegaU,OmegaV)*...
                dVdt([indOsup(2),indOsub(1)],indT(dVInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% solve the equation
A = [s(indOsub(2)),s(indOsub(1));s(indOsub(1)),s(indOsub(2))];
y = [DUikVjl-DslOmegaUkl-DskOmegaVkl;
    -DUilVjk-DskOmegaUkl-DslOmegaVkl];

DOmega = linsolve(A,y);
output = DOmega(1);

end


function [ output ] = dOmegaVdt( indOsup, indOsub, indT, U, V, s, OmegaU, OmegaV )

N = size(indT,1);

% DslOmegaUkl
DslOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaUkl = DslOmegaUkl + dsdt(indOsub(2),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaU(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DslOmegaUkl = DslOmegaUkl + dsdt(indOsub(2),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaUdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DskOmegaVkl
DskOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaVkl = DskOmegaVkl + dsdt(indOsub(1),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaV(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DskOmegaVkl = DskOmegaVkl + dsdt(indOsub(1),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaVdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

%DskOmegaUkl
DskOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaUkl = DskOmegaUkl + dsdt(indOsub(1),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaU(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DskOmegaUkl = DskOmegaUkl + dsdt(indOsub(1),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaUdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DslOmegaVkl
DslOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaVkl = DslOmegaVkl + dsdt(indOsub(2),indT,U,V,s,OmegaU,OmegaV)...
                *OmegaV(indOsup(1),indOsup(2),indOsub(1),indOsub(2));
        else
            DslOmegaVkl = DslOmegaVkl + dsdt(indOsub(2),indT(dsInd,:),U,V,s,OmegaU,OmegaV)...
                *dOmegaVdt(indOsup,indOsub,indT(dOmegaInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DUikVjl
DUikVjl = U(indOsup(1),indOsub(1))*...
            dVdt([indOsup(2),indOsub(2)],indT,U,V,s,OmegaU,OmegaV);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUikVjl = DUikVjl + dUdt([indOsup(1),indOsub(1)],indT,U,V,s,OmegaU,OmegaV)*...
                V(indOsup(2),indOsub(2));
        else
            DUikVjl = DUikVjl + dUdt([indOsup(1),indOsub(1)],indT(dUInd,:),U,V,s,OmegaU,OmegaV)*...
                dVdt([indOsup(2),indOsub(2)],indT(dVInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% DUilVjk
DUilVjk = U(indOsup(1),indOsub(2))*...
            dVdt([indOsup(2),indOsub(1)],indT,U,V,s,OmegaU,OmegaV);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUilVjk = DUilVjk + dUdt([indOsup(1),indOsub(2)],indT,U,V,s,OmegaU,OmegaV)*...
                V(indOsup(2),indOsub(1));
        else
            DUilVjk = DUilVjk + dUdt([indOsup(1),indOsub(2)],indT(dUInd,:),U,V,s,OmegaU,OmegaV)*...
                dVdt([indOsup(2),indOsub(1)],indT(dVInd,:),U,V,s,OmegaU,OmegaV);
        end
    end
end

% solve the equation
A = [s(indOsub(2)),s(indOsub(1));s(indOsub(1)),s(indOsub(2))];
y = [DUikVjl-DslOmegaUkl-DskOmegaVkl;
    -DUilVjk-DskOmegaUkl-DslOmegaVkl];

DOmega = linsolve(A,y);
output = DOmega(2);

end

