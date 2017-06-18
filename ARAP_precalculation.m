function [ mesh ] = ARAP_precalculation( mesh )
%ARAP_precalculation 
%   This function will calculate the fixed part that will be used to 
%   calculate ARAP energy and its related gradient to accelerate the 
%   moderately slow program. The results will be stored in structured 
%   input mesh variable.

%   The program organizes as follows:
%   First, calculate the aligned projection for each triangle facets
%   Second, calculate the cotangent variant for each half-edges, then
%   stored in the sparse array.

%   The following property will be added by this program:
%   mesh.vproj: Mby3by2 to project each triangle face onto 2d plane
%   respectively
%   mesh.hecot: calculate the cotangent of the angle corresponding with the
%   specific half-edge



vlen = size(mesh.v,1);
flen = size(mesh.f,1);

mesh.vproj = zeros(flen,3,2);
mesh.fcot= zeros(flen,3);
%mesh.hecot = mesh.te;
%mesh.heproj= mesh.te;

for i=1:flen
    vertid= mesh.f(i,:);
    v1= mesh.v(vertid(1),:);
    v2= mesh.v(vertid(2),:);
    v3= mesh.v(vertid(3),:);
    a=norm(v1-v2);
    b=norm(v2-v3);
    c=norm(v3-v1);
    angle=acos((a^2+c^2-b^2)/(2*a*c));
    
    %   Project onto 2d plane with rigid form.
    p1=[0,0];
    p2=[a,0];
    p3=[c*cos(angle), c*sin(angle)];
    mesh.vproj(i,1,:)= p2-p1;
    mesh.vproj(i,2,:)= p3-p2;
    mesh.vproj(i,3,:)= p1-p3;
    
    % get the cotangent value
    xcot = [vcot(v1-v3,v2-v3),vcot(v2-v1,v3-v1),vcot(v3-v2,v1-v2)];
    mesh.fcot(i,:)= xcot;
%     if mesh.te(vertid(1),vertid(2))== i
%         mesh.heproj(vertid(1),vertid(2))= p2-p1;
%     else
%         mesh.heproj(vertid(2),vertid(1))= p2-p1;
%     end
%     if mesh.te(vertid(2),vertid(3))== i
%         mesh.heproj(vertid(2),vertid(3))= p3-p2;
%     else
%         mesh.heproj(vertid(3),vertid(2))= p3-p2;
%     end
%     if mesh.te(vertid(3),vertid(1))== i
%         mesh.heproj(vertid(3),vertid(1))= p1-p3;
%     else
%         mesh.heproj(vertid(1),vertid(3))= p1-p3;
%     end
end

% calculate face area of the meshes for further use
mesh.farea= faceArea(mesh);


end

