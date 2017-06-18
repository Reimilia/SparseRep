function val = ARAP_energy_singleface( mesh, u, findex )
%ARAP_ENERGY_SINGLEFACE 此处显示有关此函数的摘要
%   此处显示详细说明

A=mesh.farea(findex);
vertid= mesh.f(findex,:);
x1= mesh.vproj(findex,1,:);
x2= mesh.vproj(findex,2,:);
x3= mesh.vproj(findex,3,:);

x1=reshape(x1,1,2);
x2=reshape(x2,1,2);
x3=reshape(x3,1,2);

xcot = mesh.fcot(findex,:);
u1= u(vertid(1),:)-u(vertid(2),:);
u2= u(vertid(2),:)-u(vertid(3),:);
u3= u(vertid(3),:)-u(vertid(1),:);
% Calculate covariance matrix
S= xcot(1)*x1'*u1 + xcot(2)*x2'*u2+ xcot(3)*x3'*u3;
% Do svd decomposition
[U,sigma,V]= svd(S);
L=U*V';
% calculate the energy
%val= A*((sigma(1,1)-1)^2+(sigma(2,2)-1)^2);
val=    xcot(1)*norm(u1-x1*L')^2+...
    xcot(2)*norm(u2-x2*L')^2+xcot(3)*norm(u3-x3*L')^2;
val=val*0.5;

end

