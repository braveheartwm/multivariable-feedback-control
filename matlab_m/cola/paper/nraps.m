function x=nraps(fil,x0);

N=max(size(x0));
p0=zeros(size(x0));
x=x0;
fold=feval(fil,x0);
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
