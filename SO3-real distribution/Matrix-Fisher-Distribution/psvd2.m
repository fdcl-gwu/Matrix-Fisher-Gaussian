function [ U, S, V ]=psvd2( F )

[U,S,V]=svd(F);
S=S*det(U)*det(V);
U=U*det(U);
V=V*det(V);

end


