% Example 8.10 RS of spinning satellite 
%    
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl8_10.m,v 1.2 2004/04/14 13:16:22 vidaral Exp $

clear all;close all;
%% Plant (3.88):
alpha = 10;
a = [0 alpha ; -alpha 0];
b = [1 0; 0 1];
c = [1 alpha; -alpha 1];
d = [ 0 0;  0 0];
G = ss(a,b,c,d);
K=eye(2);
cls=loopsens(G*K,eye(2));
w=logspace(-2,2,41);
T=frd(cls.To,w);

%% computing mu
disp('Full-block uncertainty:');
uncblk=[2 2];
[bounds,muinfo] = mussv(T,uncblk);
pk=norm(bounds(:,1),inf) % Upper bound

% Figure 8.12
bodemag(bounds(1),bounds(2)) %both upper and lower bound
xlabel('Frequency'); ylabel('Magnitude');
axis ([0.01 100, 1e-2,20]); hold on; drawnow;
plot([0.01 100],[1 1],'k:')

disp('Diagonal uncertainty:');
uncblk=[1 1; 1 1];
[bounds,muinfo] = mussv(T,uncblk);
xlabel('Frequency'); ylabel('Magnitude');
axis ([0.01 100, 1e-2,20]);hold on;drawnow;
plot([0.01 100],[1 1],'k:')
pk=norm(bounds(:,1),inf) % Upper bound

disp('Repeated complex scalar uncertainty:');
uncblk=[2 0];
[bounds,muinfo] = mussv(T,uncblk);
xlabel('Frequency'); ylabel('Magnitude');
axis ([0.01 100, 1e-2,20]);hold on;drawnow;
plot([0.01 100],[1 1],'k:')
pk=norm(bounds(:,1),inf) % Upper bound

disp('Repeated real scalar uncertainty:');
uncblk=[-2 0];
[bounds,muinfo] = mussv(T,uncblk);
xlabel('Frequency'); ylabel('Magnitude');
axis ([0.01 100, 1e-2,20]);hold on;drawnow;
plot([0.01 100],[1 1],'k:')
pk=norm(bounds(:,1),inf) % Upper bound
