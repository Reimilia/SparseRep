function [ val ] = ARAP_energy( mesh, u )
%ARAP_ENERGY 
% Given the input mesh and its parameterization
% We calculate the ARAP energy
%
%%
% $\sum_t A_t((\sigma_{1,t}-1)^2+(\sigma_{2,t}-1)^2)$

flen= length(mesh.f);

% Now we use the covariance to replace jacobian and get $\sigma$ by using
% SVD decomposition for it

val=0;

for i=1:flen
    val=val+ARAP_energy_singleface(mesh,u,i);
end

%val=val*0.5;

end

