
% PI control of unstable plant

clear all
% G(s) = 2/(s-1) (0.02s+1)^2
%
s=tf('s');
% den2 = conv([0.02 1],[0.02 1]);
% deng=conv([1 -1],den2);
% numg=4;
% G=nd2sys(numg,deng);   %Plant
G=4/((s-1)*(0.02*s+1)^2);

Kc=1.25;taui=1.5;
% numc=Kc*conv([1 1/taui],[taud 1]); 
% denc=conv([1 0],[0.000001*taud 1]);
% K=nd2sys(numc,denc); 
K=Kc*(1+1/(taui*s));
 % PI-controller (which is almost identical to 
		% the S/KS H-infinity controller obtained with the weight in Table 2.3:
        % M=1.5; wb=10; Wp = nd2sys([1/M wb], [1 1.e-6]); Wu=0.2; ) 

% L = mmult(G,K);
% s = minv(madd(1,L));
% t = msub(1,s);
L=G*K;
S=inv(1+L);
T=1-S;

%graphics;
% Step response
% [A,B,C,D]=unpck(t);[y1,x,ty]=step(A,B,C,D);
[y1,ty]=step(T);
%Input signal
% [A,B,C,D]=unpck(mmult(K,s)); [u1,x]=step(A,B,C,D,1,ty);
Us=K*S;
[u1]=step(Us,ty);
40000000000000000
figure(1)
33300000000000000
subplot(2,1,1);plot(ty,y1,ty,u1,'--',[0 4],[1 1],':' ) 
20000000000000
text(2,1.3,'y(t)')
text(2,-0.1,'u(t)')
xlabel('Time [sec]')
% Note that there is an inverse response in the input for an unstable plant

% Bode plots of L, S, T
omega=logspace(-1,2,201);
[Lm,Lp] = bode(L,omega);
Lm2(1:length(Lm))=Lm(1,1,:); Lp2(1:length(Lp))=Lp(1,1,:);
[Sm,Sp] = bode(S,omega);
Sm2(1:length(Sm))=Sm(1,1,:); Sp2(1:length(Sp))=Sp(1,1,:);
[Tm,Tp] = bode(T,omega);
Tm2(1:length(Tm))=Tm(1,1,:); Tp2(1:length(Tp))=Tp(1,1,:);


200000000000000000
figure(2)
2000000000
subplot(2,1,1);loglog(omega,[Lm2; Sm2; Tm2 ],omega,ones(size(Lm2)),':')
200000
axis([.1,100,.05,10])
text(10,1.5,'|S|'); text(2.3,3,'|L|'); text(0.5,1.5,'|T|')
%vplot('bode_l',[.1 1000,.1,5],[.1,1000,-300,300],L_f,s_f,t_f,1,'w:')
% bode(L,S,T,tf(1,1),':',omega)
% figure(2);clf;subplot(2,1,1);
% vplot('liv,lm',L_f,s_f,t_f,1,'w:'); axis([.1,100,.05,10])
% semilogx(omega,L_f)
% text(10,1.5,'|S|'); text(2.3,3,'|L|'); text(0.5,1.5,'|T|')
subplot(2,1,2);
% pl=180/pi*unwrap(angle(L_f));pl=sel(pl,1:length(omega),1)-360;
% ps=180/pi*unwrap(angle(s_f));ps=sel(ps,1:length(omega),1);
% pt=180/pi*unwrap(angle(t_f));pt=sel(pt,1:length(omega),1);
% semilogx(omega,pl,omega,ps,omega,pt,omega,0,'w:',omega,-180,'w:');
semilogx(omega,[Lp2-360; Sp2; Tp2],omega,0,':',omega,-180,'w')
axis([.1,100,-400,100])
text(.25,-290,'phaseL');
text(.25,-90,'phaseS');
text(.25,40,'phaseT');
xlabel('Frequency')

% Find GM, PM and the peak values of $s$ and $t$ 
[A,B,C,D]=ssdata(L); [Gm,Pm,W180,Wc]=margin(A,B,C,D)
norm(S,inf)  % no peak in S
norm(T,inf)  % peak in T of 1.4 at frequency 1.26 rad/s

% Summary:
% For unstable plants the peak in T is usually larger and at a lower
% frequency than that for S, and this is confirmed in this example. 
% Note that the opposite is usually true for a stable plant.  
% deadtime=0;
% [Gm,Pm,Ms,Mt]=marginsminbode(numg,deng,deadtime,Kc,taui,taud)
[Gm,Pm,Ms,Mt]=margin(L)

