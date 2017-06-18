function [ grad ] = ARAP_energy_gradient( mesh, u )
%ARAP_ENERGY_GRADIENT 
% The function is used to get the numerical gradient for ARAP
% energy based on direct calculation.
% The gradient is related with 1-ring area for each vertices

vlen = size(mesh.v,1);
flen= size(mesh.f,1);
grad= zeros(vlen,2);

% Calculate Auxiliary variable , i.e. the affine transformation part in
% the equation
for i=1:flen
    vertid= mesh.f(i,:);
    x1= mesh.vproj(i,1,:);
    x2= mesh.vproj(i,2,:);
    x3= mesh.vproj(i,3,:);
    
    x1=reshape(x1,1,2);
    x2=reshape(x2,1,2);
    x3=reshape(x3,1,2);
    
    xcot = mesh.fcot(i,:);
    u1= u(vertid(2),:)-u(vertid(1),:);
    u2= u(vertid(3),:)-u(vertid(2),:);
    u3= u(vertid(1),:)-u(vertid(3),:);
    % Calculate covariance matrix
    S= xcot(1)*u1'*x1 + xcot(2)*u2'*x2+ xcot(3)*u3'*x3;
    % Do svd decomposition
    [U,~,V]= svd(S);
    % Set the signular value all to 1
    L=U*V';
    % calculate the gradient related with the triangle
   
    grad(vertid(1),:) = grad(vertid(1),:)+ xcot(1)*(-u1) + xcot(3)*u3 +...
        xcot(1)*x1*L' + xcot(3)*(-x3)*L';
    grad(vertid(2),:) = grad(vertid(2),:)+ xcot(2)*(-u2) + xcot(1)*u1 +...
        xcot(2)*x2*L' + xcot(1)*(-x1)*L';    
    grad(vertid(3),:) = grad(vertid(3),:)+ xcot(3)*(-u3) + xcot(2)*u2 +...
        xcot(3)*x3*L' + xcot(2)*(-x2)*L';
    
end

% 
% for i=1:len
%     Get the vertice id for one ring faces
%     vring_id = oneringv(mesh);
%     k= size(vring_id,1);
%     for j=1:k
%         vid= vring_id(j);
%         face_id1= mesh.te(i,vid);
%         face_id2= mesh.te(vid,i);
%         grad(:,i)=grad(:,i)+(mesh.hecot(i,vid)+mesh.hecot(vid,i))*...
%             (u(i,:)-u(vid,:))+(mesh.hecot(i,vid)*L(face_id1,:,:)+...
%             mesh.hecot(vid,i)*L(face_id2,:,:))*(mesh.heproj(i,vid)-mesh.hevproj(vid,i));
%     end
% end
grad=reshape(grad,1,2*vlen);

end

