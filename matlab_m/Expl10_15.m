



% Example 10.15
s=tf('s');
G=[1 0; 5 1];
K=1/s*[1 0; -5 1];
S=inv(eye(2)+G*K);
T=G*K*S;

t=0:.01:40;
r1=ones(1,length(t));
r2=[zeros(1,(length(t)-1)/2) ones(1,(length(t)+1)/2)];
r=[r1;r2];
y=lsim(T,r,t);
figure
plot(t,y,t,r,':')
xlabel('time'); ylabel('y(t)');
axis([0 40 -.5 1.5])
text(4,.9,'y_1'); text(22,.7,'y_2'); text(20,1.1,'setpoint')


%Adding uncertainty
Gu=[1.2 0; 0 0.8];
Gn=Gu*G;
Sn=inv(eye(2)+Gn*K);
Tn=Gn*K*S;
yn=lsim(Tn,r,t);
figure
plot(t,yn,t,r,:)
xlabel('time'); ylabel('y(t)');
axis([0 40 -.5 1.5])
text(2,.9,'y_1'); text(25,.7,'y_2'); text(20,1.1,'setpoint')

figure
subplot(2,1,1)
plot(t,y,t,r,':')
xlabel('(a) Closed-loop response without uncertainty'); ylabel('y(t)');
axis([0 40 -.5 1.5])
text(4,.9,'y_1'); text(22,.7,'y_2'); text(20,1.1,'setpoint')

subplot(2,1,2)
plot(t,yn,t,r,:)
xlabel('(b) Closed-loop response with uncertainty'); ylabel('y(t)');
axis([0 40 -.5 1.5])
text(2,.9,'y_1'); text(25,.7,'y_2'); text(20,1.1,'setpoint')
