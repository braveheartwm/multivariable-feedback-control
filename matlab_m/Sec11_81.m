% Section 11.8.1      Reduction of a gas turbine aero-engine model
% 
% This script calculates the results and the plots shown in section 11.8.1
% in the book.
%
% Dependencies: file aero0.mat (aero-engine models)
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite

% $Id: Sec11_81.m,v 1.8 2004/04/15 13:33:56 vidaral Exp $

clear all;close all;

set(cstprefs.tbxprefs,'MagnitudeUnits','abs','MagnitudeScale','log') ;

% Load the engine models

load aero0
% SiS 13 Nov. 2010: Changed aero0 from old to newer Matlab state space format using the following commands:
%[a,b,c,d]=unpck(G_act);G_act=ss(a,b,c,d);
%[a,b,c,d]=unpck(G_eng);G_eng=ss(a,b,c,d); clear a b c d;
% To get original old (mu toolbox) format use: load aero0old

% System description

fprintf( '\nThis system has %i states, %i inputs, %i outputs.\n\n', ...
            size( G_eng.b, 1 ), size( G_eng.b, 1 ), size( G_eng.c, 1 ) ) ;

				
G=G_eng;

[Gb,hsig]=balreal(G); % Balanced realisation

fprintf( '\nHankel singular values in Table 11.1:\n\n' ) ;

for i=1:size(hsig,1)
    fprintf('%i) %i\n', i, hsig(i)) ;
end


% Figure 11.1
% a. Model reduction to 6 states using Rezidualization.
% Option 'MatchDC' in function modred is the default value,
% and can therefore be left out.(残化）

Gr=modred(Gb,7:size(Gb.a,1),'MatchDC');  % Eliminate all states from 7th
figure(1); subplot(2,3,1);
sigma(G,'k-',Gr,'r--',{1e-2,1e4}); grid;
title('Residualization') ;
axis([1e-2 1e4 1e-6 1e2]);


% b. Model reduction to 6 states using truncation.（截项）

Gt=modred(Gb,7:size(Gb.a,1),'truncate');
subplot(2,3,2);
sigma(G,'k-',Gt,'r--',{1e-2,1e4}); grid ;
title('Truncation') ;
axis([1e-2 1e4 1e-6 1e2]);


% c. Model reduction using Optimal hankel norm approximation.

[Gh,HankInfo]=hankelmr(Gb,6);

% Alternative: use ohklmr.m
%[Gh,BoundHan,Hsingval]=ohklmr(Gb,1,6);

subplot(2,3,3);
sigma(G,'k-',Gh,'r--',{1e-2,1e4}); grid ;
title('Hankel') ;
axis([1e-2 1e4 1e-6 1e2]);
% return
% Generate error systems.

Er=G-Gr;
Et=G-Gt;
Eh=G-Gh;


% Figure 11.2

figure(2); subplot(2,2,1);
sigma(Er,':',Et,'r-',Eh,'g--'); grid;
axis([1e-2 1e4 1e-10 10]);

[resE,resEf]=norm(Er,inf);
[truE,truEf]=norm(Et,inf);
[HanE,HanEf]=norm(Eh,inf);

fprintf( '\nThe maximum error for residualization is %0.4f, and occurs at %0.1f rad/s.\n', resE, resEf ) ;
fprintf( 'The maximum error for truncation is %0.4f, and occurs at %0.1f rad/s.\n', truE, truEf ) ;
fprintf( 'The maximum error for Hankel-norm optimal approximation is %0.4f, and occurs at %0.1f rad/s.\n', HanE, HanEf ) ;

BoundResTru=2*sum(hsig(7:15));

fprintf( '\nThe upper bound for residualization and truncation (\"twice the sum of the tail\") is %0.4f;\n', BoundResTru ) ;
fprintf( 'the upper bound for Hankel-norm optimal approximation is:' ) ; HankInfo


% Scaled truncation.

G0 =freqresp(G,0);
Gt0=freqresp(Gt,0);
Ws =inv(Gt0)*G0;
Gts=Gt*Ws;


% Scaled optimal Hankel-norm approximation.

Gh0=freqresp(Gh,0);
Ghs=Gh*Ws;


% Error systems with scalings.

Ets=G-Gts;
Ehs=G-Ghs;

subplot(2,2,2);
sigma(Er,':',Ets,'r-',Ehs,'g--'); grid ;
axis([1e-2 1e4 1e-10 10]);

[struE,struEf]=norm(Ets,inf);
[sHanE,sHanEf]=norm(Ehs,inf);
fprintf( 'The maximum error for scaled truncation is %0.4f, and occurs at %0.1f rad/s.\n', struE, struEf ) ;
fprintf( 'The maximum error for scaled Hankel-norm optimal approximation is %0.4f, and occurs at %0.1f rad/s.\n', sHanE, sHanEf ) ;


% Impulse responses, Figure 11.3

% Real system:

t=0:0.0001:0.03;
[y,t,x]=impulse(G,t);

% Balanced residualization

[yr,t,x]=impulse(Gr,t);

figure(3); subplot(2,3,1)
set(gcf,'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1])
plot(t,y(:,:,2),'-',t,yr(:,:,2),'--'); 
axis([ 0 0.03 -190 150 ])
xlabel('Time [s]');

% Scaled Balanced Truncation

[yt,t,x]=impulse(Gts,t);

subplot(2,3,2);
set(gcf,'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1])
plot(t,y(:,:,2),'-',t,yt(:,:,2),'--');
axis([0 0.03 -190 150])
xlabel('Time [s]');

% Scaled Hankel-norm approximation

[yh,t,x]=impulse(Ghs,t);

subplot(2,3,3);
set(gcf,'DefaultAxesColorOrder', [ 1 0 0; 0 1 0; 0 0 1 ] )
plot(t,y(:,:,2),'-',t,yh(:,:,2 ),'--');
axis([ 0 0.03 -190 150 ])
xlabel('Time [s]');


% Step responses, Figure 11.4

% Real system

t=0:0.005:1;
[y,t,x]=step(G,t);

% Balanced residualization

[yr,t,x]=step(Gr,t);

figure(4); subplot(2,3,1)
set(gcf,'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1])
plot(t,y(:,:,2),'-',t,yr(:,:,2),'--'); 
axis([ 0 1 -17.5 2.5 ]);
xlabel('Time [s]');

% Scaled Balanced Truncation

[yt,t,x]=step(Gts,t);

subplot(2,3,2);
set(gcf,'DefaultAxesColorOrder', [ 1 0 0; 0 1 0; 0 0 1 ] )
plot(t,y(:,:,2),'-',t,yt(:,:,2),'--'); 
axis([ 0 1 -17.5 2.5 ]);
xlabel('Time [s]');

% Scaled Hankel-norm approximation

[yh,t,x]=step(Ghs,t);

subplot(2,3,3);
set(gcf,'DefaultAxesColorOrder', [ 1 0 0; 0 1 0; 0 0 1 ] )
plot(t,y(:,:,2),'-',t,yh(:,:,2 ),'--'); 
axis([ 0 1 -17.5 2.5]);
xlabel('Time [s]');
