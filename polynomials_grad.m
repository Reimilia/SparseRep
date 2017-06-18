function val = polynomials_grad(c,p_basis,x)
    % convert the basis function to the whole polynomials and get its gradient
    % value at the point x
global p_basis_dx
global p_basis_dy
index= find(c~=0);
n= length(index);
val=[0,0];

for k=1:n
    val(1)=val(1)+c(index(k))*feval(p_basis_dx{index(k)},x(1),x(2));
    val(2)=val(2)+c(index(k))*feval(p_basis_dy{index(k)},x(1),x(2));
end

end