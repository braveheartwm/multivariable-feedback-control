function [gf,cn]=uncreg(w)
%gf=uncreg(w) produces an uncertainty region of plant (7.17)
%             at frequency w. The region gf is a complex vector
%             and can be plotted by: plot(real(gf),imag(gf)).
%
n=11;
p=linspace(2,3,n);
tau=2;theta=2;
gf=[];
for k=2:3,
gf=[gf;k*exp(-j*w*theta)/(j*w*tau+1)];
end
for l=2:n,
theta=p(l);
gf=[gf;k*exp(-j*w*theta)/(j*w*tau+1)];
end
for l=2:n,
tau=p(l);
gf=[gf;k*exp(-j*w*theta)/(j*w*tau+1)];
end
k=2;
gf=[gf;k*exp(-j*w*theta)/(j*w*tau+1)];
for l=n-1:-1:1,
theta=p(l);
gf=[gf;k*exp(-j*w*theta)/(j*w*tau+1)];
end
for l=n-1:-1:1,
tau=p(l);
gf=[gf;k*exp(-j*w*theta)/(j*w*tau+1)];
end
