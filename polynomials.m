function val = polynomials(c,p_basis,x)
    % convert the basis function to the whole polynomials and get its value
    % at the point x
index= find(c~=0);
n= length(index);
val=0;
for k=1:n
    val=val+c(index(k))*feval(p_basis{index(k)},x(1),x(2));
end

end