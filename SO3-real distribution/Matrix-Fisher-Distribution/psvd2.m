function [ U, S, V ]=psvd2( F )

[U,S,V]=svd(F);
for i = 1:3
    temp = nonzero(U(:,i));
    if temp(1) < 0
        U(:,i) = -U(:,i);
        V(:,i) = -V(:,i);
    end
end

S=S*det(U)*det(V);
U=U*det(U);
V=V*det(V);

end


