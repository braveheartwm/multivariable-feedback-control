function [t,x]= exeu( fn, tv, x0)
%
%  Function to simulate xdot = fn(t,x); using explicite euler.
%
   t0 = tv(1);   tf = tv(2);  dt = tv(3);
   
   N = ceil((tf-t0)/dt);
   nx = max(size(x0));
   
   x = zeros(N+1,nx);
   t = zeros(N+1,1);
   x(1,:) = x0.';
   t(1)   = t0;
   for i=2:N+1
      t(i) = t(i-1)+dt;
      x(i,:)= x(i-1,:)+dt*feval(fn, t(i), x(i-1,:).' ).';
   end
end    
