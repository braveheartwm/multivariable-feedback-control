clear all
x0 = [0.9; 0.4737; 0.1];
[t,x] = ode45('dist', 0,25,x0,1e-6,1);

for i = 1:size(x,1)
  dxd(i,:)=x(i,:)-x0.';
end

dx = mpck(dxd,t);
l1 = [0     4.5
      0.0168 4.5
      2     Inf];
vplot(sel(dx,1,1),'-',sel(dx,2,1),'--',sel(dx,3,1),'-.',l1,':')
xlabel('Time')
ylabel('Composition')
set(gca,'XTick',[0 4.5 10 15 20 25]);
set(gca,'FontSize',14)
text(15,0.0084,'x1','VerticalAlignment','top');
text(15,0.0261,'x2','VerticalAlignment','bottom');
text(15,0.0109,'x3','VerticalAlignment','bottom');

%print -deps st3stg.eps
%!mv st3stg.eps ../latex/st3stg.eps
