function [ tar_mesh ] = SurfaceApprox( mesh, p_basis, options )
%SURFACEAPPROX 
%   Do surface approximation by sparse representation.
%   The program run 3 times for each coordinate component
%   The parameterization process will get the initial parameter
%   For now we use simple Spectural Conformal Mapping to do this
tar_mesh=mesh;
%mesh.u = embedSCP(mesh,'fiedler');
mesh.u = mesh.v(:,1:2);

% regularization
% min_x= min(mesh.u(:,1));
% min_y= min(mesh.u(:,2));
% max_x= max(mesh.u(:,1));
% max_y= max(mesh.u(:,2));
% mid_x= (min_x+max_x)/2;
% mid_y= (min_x+max_y)/2;
% 
% Note : this is important if the result of parameterization will not be
%  embedded into [0,1] * [0,1]. Otherwise the gradient will blow up or the
%  result will suffer severe distortion :(
% mesh.u(:,1)= (mesh.u(:,1)-min_x)/(max_x-min_x);
% mesh.u(:,2)= (mesh.u(:,2)-min_y)/(max_y-min_y);

% precalculation for ARAP energy
mesh= ARAP_precalculation(mesh);

%For test purpose
figure(1);plotMesh(mesh, 'uefb');
%[ Dangle, Darea, L2, Linf, Cangle, Carea ] = meshDistortion(mesh, 1);
%printf('[SCP Fiedler] DAngle: %.4f  Darea: %.4f  L2: %.4f  Linf: %.4f  Conf: %.4f  Auth: %.4f\n', Dangle, Darea, L2, Linf, Cangle, Carea);  


approx_v= zeros(size(mesh.v,1),3);

for i=3:3
    max(abs(mesh.v(:,i)))
    [c,u]= SparseRep(mesh.v(:,i),p_basis,mesh,mesh.u,options);
    % Test plot
    norm(mesh.u-u);
    temp_mesh= mesh;
    temp_mesh.u=u;
    figure(i+1);
    plotMesh(temp_mesh,'uefb');
    for k=1:size(mesh.v,1)
        approx_v(k,i) = polynomials(c,p_basis,u(k,:));
    end
end

% Add non-linear transfer part to the generated mesh
approx_v(:,1)=u(:,1);
approx_v(:,2)=u(:,2);
tar_mesh.v= approx_v;
tar_mesh.n= estimateNormal(tar_mesh);

end

