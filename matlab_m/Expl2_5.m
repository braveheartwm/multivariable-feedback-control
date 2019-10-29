% ExplX2_5_15.m, EXTRA example, 
% PI control of unstable plant

% G(s) = 4/(s-1) (0.02s+1)^2
%G=nd2sys([4],[0.0004 0.0396 0.96 -1]);
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: ExplX2_1MOD.m,v 1.1 2004/10/13 09:38:07 kariwala Exp $

clear
G=tf([4],[0.0004 0.0396 0.96 -1]); %Plant

K = tf([1.875 1],[1.5 0]);   % PI-controller (which is almost identical to 
		% the S/KS H-infinity controller obtained with the weight in Table 2.3:
        % M=1.5; wb=10; Wp = tf([1/M wb], [1 1.e-6]); Wu=2; ) 

L = G*K;
S = 1/(1+L);
T = 1-S;

% Step response
[y1,ty]=step(T,4);
%Input signal
[u1,tu]=step(K*S,ty);
figure(1);clf;
plot(ty,y1,ty,u1,'--')   
hold
plot(ty,y1./y1,':')
xlabel('Time[s]')
axis([0 4 -0.5 1.5]);

% Bode plots of L, S, T
figure(2);clf;
[mag1,pha1,w]=bode(L);
[mag2,pha2]=bode(S,w);
[mag3,pha3]=bode(T,w);
subplot(2,1,1)
loglog(w,mag1(:),w,mag2(:),w,mag3(:))
ylabel('Magnitude')
axis([0.1 1e2 0.05 10]);
hold %on
plot(w,w./w,':')
text(0.82,1.6,'|T|'), text(1,5,'|L|'), text(0.8,0.38,'|S|') 
hold %off

subplot(2,1,2)
semilogx(w,pha1(:),w,pha2(:),w,pha3(:),w,-180*w./w,':')
axis([0.1 1e2 -400 100]);
xlabel('Frequency'),ylabel('Phase')
text(1,-35,'PhaseT'), text(1,-130,'PhaseS'), text(1,-240,'PhaseL') 

% Find GM, PM and the peak values of $s$ and $t$ 
[Gm,Pm,W180,Wc]=margin(L)
[norm_S,ws]=norm(S,inf,1e-4) % peak in S of 1.1 at frequency 21.2 rad/s
[norm_T,wt]=norm(T,inf,1e-4) % peak in T of 1.4 at frequency 1.3 rad/s
