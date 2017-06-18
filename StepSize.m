function [alfa,x] = StepSize(f,g, x, d, alfa, params)
%  Implements simple Wolfe conditions.

x0 = x.p;
Dphi0 = x.g'*d;
phi0 = x.f;
alfaL = 0;
alfaR = inf;

fprintf('initialize: Dphi0 = %f, phi0= %f\n\n',Dphi0,phi0);

if ( alfa <= 0 || Dphi0 > 0 )
  error('Initialization of step incorrect');
end;

c1 = params.ftol;
c2 = params.gtol;


iter = 0;
while abs(alfaR-alfaL) > params.xtol
  iter = iter + 1;
  x.p = x0 + alfa*d;
  x.f = feval(f,x.p);  
  fprintf('iter %d : func_val = %f, tol = %f\n ',iter,x.f,phi0+alfa*c1*Dphi0);
  
  if (x.f >= phi0 + alfa*c1*Dphi0) 
    alfaR = alfa;
    fprintf('switch 1 \n');
    alfa = Interp(alfaL, alfaR);
  else
    x.g = feval(g,x.p);
    DphiAlfa = x.g'*d;
    fprintf('switch 2: DphiAlfa=%f, tol=%f\n',DphiAlfa,c2*Dphi0);
    if (DphiAlfa >= c2*Dphi0)
        return;
    else
      alfaL = alfa;
    end
    if isinf(alfaR)
      alfa = Extrap(alfaL);
    else
      alfa = Interp( alfaL, alfaR);
    end	
  end
end
fprintf('step size criteria were not met\n');
fprintf('after %d step size iterations.\n', iter);
%error('STEP SIZE FAILURE');
end

function mid = Interp(left, right)
mid = (left + right)/2.0;
end

function mid = Extrap(left)
mid = 2*left;
end

