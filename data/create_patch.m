mesh=readMesh('plane_z.obj');

mesh.v(:,3)= 1-abs(mesh.v(:,1)+mesh.v(:,2));
mesh.n= estimateNormal(mesh);

writeMesh(mesh,'plane_bent2.obj');
plotMesh(mesh,'efb');