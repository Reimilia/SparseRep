function [ grad ] = monomial_gradients(y,c,p_basis,u )
%MONIMIAL_GRADIENTS 
%   calculate the gradient with the following formula
%
%           $\lVert y - p(x) \rVert^2$
%
%   Here p(x) is a polynomial for x

m=size(u,1);
grad=zeros(m,2);

eps=1e-4;
func= @(x) polynomials(c,p_basis,x);
for i=1:m
    % lazy way to do that
    grad(i,:)= 2*(y(i)-polynomials(c,p_basis,u(i,:)))*polynomials_grad(c,p_basis,u(i,:));
    %(func(u(i,:))-func(u(i,:)+[eps,0]))/eps;
    %grad(i,2)= polynomials_grad(c,p_basis,u(i,:));
    %(func(u(i,:))-func(u(i,:)+[0,eps]))/eps;
end
grad= reshape(grad,1,2*m);
end




