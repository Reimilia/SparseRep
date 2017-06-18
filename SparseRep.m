function [ c,U ] = SparseRep( b,p_basis,mesh,u0,options )
%SPARSEREP The main optimize solver for Sparse representation.
% Input :   b: Scalar Observations
%           p_basis: monomials (using feval in matlab to activate it), in
%           this program this shall be 2 variables: @(x,y) x^ay^b
%           mesh: input 3d mesh for this problem, use for calculation of
%           ARAP energy
%           u0: initial parameterization
%           options: selected optional settings for the optimization
% 
%
%
%   The program aims at solving C,U as proposed in paper to get the 
%   sparse approximation result

%parameter of the optimization
sparse_para=get_options(options,'s',15);
beta1 = get_options(options,'beta1',0.2);
beta2 = get_options(options,'beta2',0.1);
rho =  get_options(options,'rho',0.1);
tau =  get_options(options,'tau',0.5);

% iteration times
n_iter = get_options(options,'iterstep',20);

% parameter for debugging
test_on = get_options(options,'debug',0);
err_tolerance = get_options(options,'eps',1e-3);

% Initialize variables
lambda=zeros(1,size(u0,1));
U=u0;
c=zeros(1,length(p_basis));
b=b';
f=zeros(size(b));

loop_count=1;
err=1;
%We need to formulate the function that will be used in optimization
%process, like this:
%   func= @(x)optimization_function(x,...
%             f,c,u0,p_basis,mesh,beta1,beta2,rho,lambda);
%   func_g= @(x)optimization_function_gradient(x,...
%             f,c,u0,p_basis,mesh,beta1,beta2,rho,lambda);

while(loop_count<n_iter && err>=err_tolerance)
    P= get_numerical_basis(p_basis,U);
    % Sparsity Optimization
    last_c=c;
    c= OMP(f+lambda/rho,P,sparse_para);
    % Parameterization Optimization
    func= @(x)optimization_function(x,...
            f,c,u0,p_basis,mesh,beta1,beta2,rho,lambda);
    func_g= @(x)optimization_function_gradient(x,...
             f,c,u0,p_basis,mesh,beta1,beta2,rho,lambda);
  
    last_U=U;
    U= BFGS_gradient_approximate(reshape(U,[],1),func,func_g, []);
    U= reshape(U,[],2);
    %update f
    P= get_numerical_basis(p_basis,U);
    a= c*P'-lambda/rho;
    f= (2*b+rho*a)/(2+rho);
    %update lambda
    lambda= lambda+ tau*rho*(f-c*P');
    %update residual error
    err= max([norm(c-last_c)/norm(c),norm(U-last_U,'fro')/norm(U,'fro'),norm(f-c*P')]);
    
    
    if test_on==1
        % do debugging or intermediate test
        fprintf('\n\nloop : %d , err = %f, c_norm = %f, f_norm = %f, approximate_err= %f \n\n', loop_count,err,norm(c),norm(f),norm(f-c*P'));
        pause(1);
    end

    loop_count= loop_count+1;
end
end

function P = get_numerical_basis(p_basis,U) 
%   U must be Mby2 matrix
    n= length(p_basis);
    m= length(U);
    P=zeros(m,n);
    
    for i=1:m
        for j=1:n
            P(i,j)= feval(p_basis{j},U(i,1),U(i,2));
            %fprintf('%f %f %f\n',P(i,j),U(i,1),U(i,2));
        end
    end
end

function [ val ] = optimization_function(u,...
        f,c,u0,p_basis,mesh,beta1,beta2,rho,lambda)
    % This function aims at wrapping the optimization function into the
    % lambda form: given the input point u
    %
    %
    u=reshape(u,[],2);
    PP= get_numerical_basis(p_basis,u);
    approx_norm= 0.5*rho*norm(f-c*PP'+lambda/rho)^2;
    Q1U= beta1*(norm(u-u0,'fro')^2);
    Q2U= beta2*ARAP_energy(mesh,u);
    fprintf('Energy Contribution: approx_norm = %f Q1U = %f Q2U = %f \n',approx_norm,Q1U,Q2U);
    %pause(2);
    val= approx_norm+Q1U+Q2U;
end

function [ val ] = optimization_function_gradient(u,...
        f,c,u0,p_basis,mesh,beta1,beta2,rho,lambda)
%   This function aims at getting the gradient of the target function
%   Then to use Quasi-Newton method with BFGS adjustment

%   Here we use the simplest difference method to approximate numerical
%   value for the gradient since this function is hard to calculate the
%   explicit formula (although it is possible to write it down)

%   Input:  u is the input point
%           p_basis is the basis function
%           mesh is the input mesh for this problem, it is mainly used to
%           calculate ARAP energy.
    u=reshape(u,[],2);
    approx_grad= 0.5*rho*monomial_gradients(f+lambda/rho,c,p_basis,u);
    Q1_grad=beta1*2*reshape(u-u0,1,2*size(u,1));
    Q2_grad=beta2*ARAP_energy_gradient(mesh,u);
    val=approx_grad + Q1_grad + Q2_grad;
    fprintf('Gradient Contribution: approx_norm = %f Q1 = %f Q2 = %f \n',...
        norm(approx_grad),norm(Q1_grad),norm(Q2_grad));
    %pause(2);
    val=val';
end