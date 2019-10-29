s = tf('s');
G1 = (-0.6*s + 1)/(6*s + 1)*tf(1,1,'iodelay',1);
G2 = 1/(6*s + 1)/(0.4*s + 1);

% Without cascade
K = 0.9*(1 + 1/(9*s));
CL = G1*G2*K*inv(1 + pade(G1*G2*K,10));
CLGd = inv(1 + pade(G1*G2*K,10))*G1;

% With cascade
K1 = 1.36*(1 + 1/(6*s));
K2 = 10.33*(1 + 1/(2.4*s));
CLcas = G1*G2*K1*K2*inv(1 + G2*K2*(1 + pade(G1*K1,10)));
CLcasGd = G1*inv(1 + G2*K2*(1 + pade(G1*K1,10)));

t = 0:0.005:100;
r(1:length(t)) = 1;
d(1:round(length(t)/2)) = 0;
d(round(length(t)/2)+1:length(t)) = 6;

y = lsim([CL CLGd],[r' d']',t);
ycas = lsim([CLcas CLcasGd],[r' d']',t);

plot(t,y,'b:',t,ycas,'r',t,r,'k-')
axis([0 100 -0.2 5])
text(85,1.3,'Setpoint'),text(28,3,'Without Cascade'),text(54,1.5,'With Cascade')
text(20,0.5,'Setpoint Change'), text(60,0.5,'Disturbance Change')
