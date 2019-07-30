function [ U, S, V ] = usvd( F, positive )

if ~exist('positive','var')
    positive = false;
end

[U,S,V]=svd(F);

for i = 1:3
    [~,ind] = max(abs(V(:,i)));
    if V(ind,i) < 0
        V(:,i) = -V(:,i);
        U(:,i) = -U(:,i);
    end
end

if positive
    U=U*diag([1 1 det(V)]);
    V=V*diag([1 1 det(V)]);
else
    S=S*diag([1 1 det(U*V)]);
    U=U*diag([1 1 det(U)]);
    V=V*diag([1 1 det(V)]);
end

end

