clear all

load TESYS

%Select the subsystem
outp = [7,8,9,11,12,13,15,16,18,21,22]; 
Gw   = Gusc(outp,1:12);

% Scale inputs by 10%; 
Du = 10*eye(12);

%Scaling for Outputs; 
De = inv(diag([54.1 1.5 1.2 1 1 52.6 1 62 1 0.2 0.2])); 
Gs = De*Gw*Du;

p1 = pole(Gs);
p = sort(p1(find(real(p1) >= 0)));

nx = size(Gs.A,1);

for i = 1:length(p)
    [U,S,V] = svd(Gs.A-p(i)*eye(nx));
    Yp(:,i) = Gs.C*V(:,nx);
    Up(:,i) = (Gs.B)'*U(:,nx);
end


