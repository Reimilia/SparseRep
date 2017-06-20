function [ u, val ] = BFGS_gradient_approximate( u0, f,g, options )
%BFGS_gradient_approximate 
%   This is the Quasi-Newton method with BFGS approximation for Hesse
%   matrix, it receives input function and its gradient, and begin
%   searching at the initial point u0. f will be the original function
%   while g will be the target function.

max_iter= get_options(options,'iterstep',10);
tolerance = get_options(options,'tol', 1e-4);

%  Initialize parameter structure for StepSize function call.
params = struct('ftol', 0.3, 'gtol', 0.4, 'xtol', 1e-4, 'stpmin', 0, ...
                'stpmax', 20, 'maxfev', 10000);

xc.p= u0;
            
I = eye(size(xc.p, 1));  % Locally stored identity matrix.
for i = 1:max_iter
    %  Compute function and gradient at current point.
    xc.f = feval(f, xc.p);
    xc.g = feval(g, xc.p);
    %  Check for termination condition: (scaled) norm of gradient less
    %  than toler.
    if norm(xc.g) < tolerance
%         inform.status = 1;  % Indicates success.
%         inform.iter = i;  % Number of iterations.
        u= xc.p;
        val = xc.f;
        return;
    end
    
    %  For the first step, we use the identity matrix as an initial inverse
    %  Hessian approximation.
    if i == 1
        H = I;
    else
        %  Update the current inverse Hessian approximation, H.
        H = (I - rho*s*y') * H * (I - rho*y*s') + rho*(s*s');
    end
    
    %  Compute the current search direction.
    
    p = -H * xc.g;
    
    %  Get step size that satisfies simple Wolfe conditions.
    %  NOTE:  alfa = 1 should always be tried first since this step length will
    %         eventually always be accepted (under certain conditions), thereby
    %         producing superlinear convergence of the overall algorithm.
    [alfa, ~] = StepSize(f,g, xc, p, 1, params);
    
    % In this case the step search failed, at some scenario, it might be 
    % the case that we find the Stationary Point
    if alfa<params.xtol
         break
    end
    %  Update current point in p-direction with step size alpha.
    xc.p = xc.p + alfa * p;

    %  Update parameters.
    s = alfa * p;  %  s = new_point - prev_point = alfa * p
    y = feval(g,xc.p) - xc.g;  % y = grad(new_point) - grad(prev_point)
    %s=s';
    %y=y';
    rho = 1 / (y'*s);
    
    if i == 1
        H = (s'*y)/(y'*y) * I;
    end
end

fprintf('\n\nWarning: BFGS exit without reaching the stop creteria!\n\n');
u=xc.p;
val= xc.f;
end

