addpath(genpath('external\matlabmesh'));

mesh=readMesh('data\saddle.obj');

% maximum degree for the monomial
degree= 14;

%plot input test mesh
%figure(1);
%plotMesh(mesh,'efb');


%Construct monomial basis
cnt=0;
basefunc=@(x,y,a,b) x^a*y^b;
p_basis={};
% lazy way to pass the gradient
% TODO: it is not a good idea to pass global variable.
global p_basis_dx;
global p_basis_dy;
p_basis_dx={};
p_basis_dy={};
for i=0:degree
    for j=0:degree-i
        cnt= cnt+1;
        if i==0&&j==0
            p_basis{cnt}=@(x,y) 1;
            p_basis_dx{cnt}= @(x,y) 0;
            p_basis_dy{cnt}= @(x,y) 0;
        else    
            p_basis{cnt}= @(x,y) basefunc(x,y,i,j);
            if i==0
                p_basis_dx{cnt}=@(x,y) 0;
            else
                p_basis_dx{cnt}=@(x,y) i*basefunc(x,y,i-1,j);
            end
            if j==0
                p_basis_dy{cnt}=@(x,y) 0;
            else
                p_basis_dy{cnt}=@(x,y) j*basefunc(x,y,i,j-1);
            end
        end
    end
end\

%Set parameters for the optimization
options=[];
% open the debug output
options.debug=1;
% change iter step here
options.iterstep=15;

% let's solve the problem
new_mesh= SurfaceApprox(mesh,p_basis,options);

% Time to plot error map
figure(10); 
plotMesh(new_mesh,'efb');
figure(11);
v_len= size(mesh.v,1);
v_err= zeros(1,v_len);
for i=1:v_len
    v_err(i)= norm(mesh.v(i,:)-new_mesh.v(i,:));
end

plotMesh(mesh,'efb',v_err');
