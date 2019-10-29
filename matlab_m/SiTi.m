s = tf('s');
k = 0.1;
z=5;p=4.5;a=10; 
g0 = [(s-z)/(s-p),k;(s-z)/(sqrt(k)*s+1),(sqrt(k)*s+1)/(s-p)];

g0 = [(s-z)/(s-p),1;0.01*(s-z)/(s+a),0.01];

X = diag([1 1]);
G = g0*X;

gp=evalfr(G,p+0.0001); [U,S,V]=svd(gp);
yp=U(:,1); up = V(:,1);

gz=evalfr(G,z); [U,S,V]=svd(gz); 
yz=U(:,length(yp)); uz = V(:,length(yp));

S = yp'*yz
SI = up'*uz


