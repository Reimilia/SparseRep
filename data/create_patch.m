mesh=readMesh('plane_z.obj');

mesh.v(:,3)= abs(mesh.v(:,1));
mesh.n= estimateNormal(mesh);

writeMesh(mesh,'plane_bent3.obj');
plotMesh(mesh,'efb');