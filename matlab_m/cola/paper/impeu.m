function [t,x] = impeu(fn,tv,x0)

   t0 = tv(1);   tf = tv(2);  dt = tv(3);
   
   N  = ceil((tf-t0)/dt);
   nx = max(size(x0));
   
   x = zeros(N+1,nx);
   t = zeros(N+1,1);
   x(1,:) = x0.';
   t(1)   = t0;
   for i=2:N+1
     
      t(i) = t(i-1)+dt;
     
      p0=zeros(size(x0));
      x=x(i,:);
      fold= feval(fn,x0);
      normold=norm(fold);
delta=0.001;



while normold>1.e-10

  for i=1:N
    p=p0; p(i)=delta;
    J(:,i)=(feval(fil,x+p)-feval(fil,x-p))/(2*delta);
  end %for

  dx=-J\fold;
  f=feval(fil,x+dx);

  while norm(f)>normold
    dx=dx/2;
    f=feval(fil,x+dx);
  end %while

  x=x+dx;
  fold=f;
  normold=norm(fold);

end %while

      
      x(i,:)= x(i-1,:)+dt*feval(fn, t(i), x(i-1,:).' ).';
      
      
      
      
   end
end    
