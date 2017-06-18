function grad= ARAP_energy_numgradient(mesh, u)


len = size(u,1);
grad=zeros(len,2);
eps=1e-4;
fun = @(x,y) ARAP_energy_singleface(mesh,x,y);
for i=1:len
    face_id= oneringf(mesh,i);
    uh=u;
    uh(i,1)=uh(i,1)+eps;
    E=0;
    for k=1:length(face_id)
        E= E+fun(uh,k)-fun(u,k);
    end
    grad(i,1)= E/eps;
    uh=u;
    uh(i,2)=uh(i,2)+eps;
    E=0;
    for k=1:length(face_id)
        E= E+fun(uh,k)-fun(u,k);
    end
    grad(i,2)= E/eps;
    
end

grad=reshape(grad,1,2*len);

end