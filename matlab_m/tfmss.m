function [A,B,C,D] = tfmss(d,h,k)
[m,r] = size(k);
I = eye(m);
n = size(d,2);
v = ones(1,(n-1)*m);
A = [];
B =[];C =[];D=[];den=[];

for i = 1:n
    di = d(1,i);
    deni = -di*I;
    den = [den;deni];
end
A = diag(v);
A1 = zeros(m,(n-1)*m);
A = [A1;A];
A = [A,den];
B = h;
C = zeros(m,(n-1)*m);
C=[C,I];
D = k;
[A,B,C,D] = minreal(A,B,C,D);
end