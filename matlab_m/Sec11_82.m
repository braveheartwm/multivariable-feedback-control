% Section 11.8.2        Reduction of an aero-engine controller
% 
% This script calculates the results and the plots shown in section 11.8.2
% in the book.
%
% Dependencies: file aeroK.mat (aero-engine and controller models)
%
% Known problems at commitment time:function sigma does not plot the 
% singular values of simple transfer functions, though it works for more
% complicated ones.
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Sec11_82.m,v 1.6 2004/02/23 14:01:58 zenith Exp $


clear all
close all
set( cstprefs.tbxprefs, 'MagnitudeUnits', 'abs', 'MagnitudeScale', 'log' ) ;


% The shaped plant of the model and the unscaled controller 
% are saved in aeroK.mat from the design, in Sec12_33.m.
% Also the reference model is inside this file.

load aeroK

% Controller K, Shaped plant Gs and reference model Mo.

%w = logspace(-4, 4, 201);

K1 = K( 1:3, 1:3 ) ;
K2 = K( 1:3, 4:6 ) ;

G = Gs;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'M' ;
cleanupsysic = 'yes' ;
sysic ;

% Figure 11.5

figure( 1 ) ; subplot( 2, 2, 1 ) ;
sigma( Mo, '-', M, '--', { 1e-2, 1e3 } ) ; grid ;  % Mo is not plotted, unknown reasons
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

[ Merr, Mpeak ] = norm( ( M - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe model error is %0.3f and occurs at %0.1f rad/s.\n', Merr, Mpeak ) ;

% Scale the prefilter to match Mo at steady state.

Mo0 = freqresp( Mo, 0 ) ;
M0  = freqresp( M, 0 ) ;

SM = inv(M0) * Mo0 ;

K1s = K1 * SM ;

systemnames  = 'G K1s K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1s + K2 ]' ;
input_to_K1s = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Ms' ;
cleanupsysic = 'yes' ;
sysic;


subplot( 2 ,2, 2 ) ;
sigma( Mo, '-',  Ms, '--', { 1e-2 1e3 } ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

fprintf( '\n\n***Results for scaled-controller approach:\n' ) ;
[ Merrs, Mpeaks ] = norm( ( Ms - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe scaled-model error is %0.2f and occurs at %0.0f rad/s.\n', Merrs, Mpeaks ) ;


% 1. Reducing the scaled controller.


% a. Balanced residualization.

Ks = [ K1s, K2 ] ;
[ Ksb Hs ] = balreal( Ks ) ;
Ksbr = modred( Ksb, 8:size( Ksb.a, 1 ) ) ;

Rbound = sum( 2 * Hs( 8:max( size( Hs ) ) ) ) ;
fprintf( '\nThe controller error bound for residualization is %0.3f.', Rbound ) ;

[ RKnorm, RKpeak  ] = norm( ( Ksbr - Ks ), inf, 1e-7 ) ;
fprintf( '\nThe actual controller error for residualization is %0.3f.\n', RKnorm ) ;

K1 = Ksbr( 1:3, 1:3 ) ;
K2 = Ksbr( 1:3, 4:6 ) ;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Mbr' ;
cleanupsysic = 'yes' ;
sysic ;

%Figure 11.6

figure( 2 ) ; subplot( 2, 3, 1 ) ;
sigma( Mo, '-', Mbr,'--', { 1e-2, 1e3 } ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

[ RMnorm, RMpeak  ] = norm( ( Mbr - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe maximum model error for residualization is %0.4f, and occurs at %0.1f rad/s.\n', RMnorm, RMpeak ) ;


% b. Balanced truncation.

Ksbt = modred( Ksb, 8:size( Ksb.a, 1 ), 'truncate' ) ;

K1 = Ksbt( 1:3, 1:3 ) ;
K2 = Ksbt( 1:3, 4:6 ) ;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Mbt' ;
cleanupsysic = 'yes' ;
sysic ;

% Scaling the truncated-controller system for steady-state invariance

Mbt0 = freqresp( Mbt, 0 ) ;
SMt  = inv( Mbt0 ) * Mo0 ;
Mbt  = Mbt * SMt ;

subplot( 2, 3, 2 ) ;
sigma( Mo, '-', Mbt, '--' ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

[ TMnorm, TMpeak  ] = norm( ( Mbt - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe maximum model error for scaled truncation is %0.2f, and occurs at %0.0f rad/s.\n', TMnorm, TMpeak ) ;


% c. Optimal hankel norm approximation.

[ Ksbh Hinfo ] = hankelmr( Ksb, 7 ) ;

fprintf( '\nThe controller error bound for unscaled optimal Hankel-norm approximation is:' ) ; Hinfo.ErrorBound

K1 = Ksbh( 1:3, 1:3 ) ;
K2 = Ksbh( 1:3, 4:6 ) ;

[ HKnorm, HKpeak  ] = norm( ( [ K1 K2 ] - Ks ), inf, 1e-7 ) ;
fprintf( 'The actual controller error for scaled optimal Hankel-norm approximation is %0.3f.\n', RKnorm ) ;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Mbh' ;
cleanupsysic = 'yes' ;
sysic ;

% Scaling the system with optimal Hankel-norm approximated controller for steady-state invariance

Mbh0 = freqresp( Mbh, 0 ) ;
SMh  = inv( Mbh0 ) * Mo0 ;
Mbh  = Mbh * SMh ;

subplot( 2, 3, 3 ) ;
sigma( Mo, '-', Mbh, '--' ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

[ HMnorm, HMpeak  ] = norm( ( Mbh - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe maximum model error for scaled optimal Hankel-norm approximation is %0.2f, and occurs at %0.0f rad/s.\n', HMnorm, HMpeak ) ;


% 2. Reduce the unscaled controller.


% a. Balance residualization.

[ Kb H ]  = balreal( K ) ;
KBr = modred( Kb, 8:size( Kb.a, 1 ) ) ;

K1 = KBr( 1:3, 1:3 ) ;
K2 = KBr( 1:3, 4:6 ) ;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Mbr2' ;
cleanupsysic = 'yes' ;
sysic ;

% Scaling the system matrix

Mbr20 = freqresp( Mbr2, 0 ) ;
SM2   = inv( Mbr20 ) * Mo0 ;
Mbr2  = Mbr2 * SM2 ;

% Figure 11.10

figure( 3 ) ; subplot( 2, 3, 1 ) ;
sigma( Mo, '-', Mbr2, '--' ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

fprintf( '\n\n***Results for unscaled-controller approach:\n' ) ;
[ RMnorm2, RMpeak2  ] = norm( ( Mbr2 - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe maximum model error for residualization is %0.2f, and occurs at %0.0f rad/s.\n', RMnorm2, RMpeak2 ) ;


% b. Balanced truncation.

KBt = modred( Kb, 8:size( Kb.a, 1 ), 'truncate' ) ;

K1 = KBt( 1:3, 1:3 ) ;
K2 = KBt( 1:3, 4:6 ) ;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Mbt2' ;
cleanupsysic = 'yes' ;
sysic ;

% Scaling the system matrix

Mbt20 = freqresp( Mbt2, 0 ) ;
SMt2  = inv( Mbt20 ) * Mo0 ;
Mbt2  = Mbt2 * SMt2 ;

subplot( 2, 3, 2 ) ;
sigma( Mo, '-', Mbt2, '--' ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

[ TMnorm2, TMpeak2  ] = norm( ( Mbt2 - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe maximum model error for scaled truncation is %0.1f, and occurs at %0.1f rad/s.\n', TMnorm2, TMpeak2 ) ;


% c. Optimal hankel norm approximation.

[ KBh Hinfo2 ]= hankelmr( Kb, 7 ) ;

fprintf( '\nThe controller error bound for unscaled optimal Hankel-norm approximation is:' ) ; Hinfo2.ErrorBound

K1 = KBh( 1:3, 1:3 ) ;
K2 = KBh( 1:3, 4:6 ) ;

systemnames  = 'G K1 K2' ;
inputvar     = '[ r( 3 ) ]' ;
outputvar    = '[ G ]' ;
input_to_G   = '[ K1 + K2 ]' ;
input_to_K1  = '[ r ]' ;
input_to_K2  = '[ G ]' ;
sysoutname   = 'Mbh2' ;
cleanupsysic = 'yes' ;
sysic ;

% Scaling the system matrix

Mbh20 = freqresp( Mbh2, 0 ) ;
SMh2  = inv( Mbh20 ) * Mo0 ;
Mbh2  = Mbh2 * SMh2 ;

subplot( 2, 3, 3 ) ;
sigma( Mo, '-', Mbh2, '--' ) ; grid ;
axis( [ 1e-2 1e3 1e-4 10 ] ) ;

[ HMnorm2, HMpeak2  ] = norm( ( Mbh2 - Mo ), inf, 1e-7 ) ;
fprintf( '\nThe maximum model error for scaled optimal Hankel-norm approximation is %0.1f, and occurs at %0.1f rad/s.\n', HMnorm2, HMpeak2 ) ;


% Step responses.


t = 0:0.01:2 ;
SimStepAxis = [ 0 2 -0.5 1.5 ] ;


% 1. Full controller.

yr1 = lsim( Ms, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
yr2 = lsim( Ms, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
yr3 = lsim( Ms, [ 0 0 1 ]' * ones( size( t ) ), t ) ;


% 2. Balanced residualization.

ybrr1 = lsim(Mbr, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
ybrr2 = lsim(Mbr, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
ybrr3 = lsim(Mbr, [ 0 0 1 ]' * ones( size( t ) ), t ) ;


% Figure 11.7

figure( 4 ) ; subplot( 2, 3, 1 ) ;
set( gcf, 'DefaultAxesColorOrder', [ 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1 ] ) ;
plot( t, yr1, '-', t, ybrr1, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 2 ) ;
plot( t, yr2, '-', t, ybrr2, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 3 ) ;
plot( t, yr3, '-', t, ybrr3, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;


% 3. Balanced truncation.

ybtr1 = lsim( Mbt, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
ybtr2 = lsim( Mbt, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
ybtr3 = lsim( Mbt, [ 0 0 1 ]' * ones( size( t ) ), t ) ;


%Figure 11.8

figure( 5 ) ; subplot( 2, 3, 1 ) ;
set( gcf, 'DefaultAxesColorOrder', [ 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1 ] ) ;
plot( t, yr1, '-', t, ybtr1, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;
			     
subplot( 2, 3, 2 ) ;
plot( t, yr2, '-', t, ybtr2, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 3 ) ;
plot( t, yr3, '-', t, ybtr3, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;


% 4. Hankel norm approximation.

ybhsr1 = lsim( Mbh, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
ybhsr2 = lsim( Mbh, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
ybhsr3 = lsim( Mbh, [ 0 0 1 ]' * ones( size( t ) ), t ) ;


% Figure 11.9

figure( 6 ) ; subplot( 2, 3, 1 ) ;
set( gcf, 'DefaultAxesColorOrder', [ 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1 ] ) ;
plot( t, yr1, '-', t, ybhsr1, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 2 ) ;
plot( t, yr2, '-', t, ybhsr2, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 3 ) ;
plot( t, yr3, '-', t, ybhsr3, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;


% 4. Unscaled model as starting point for model reduction.


% 2. Balanced residualization.

ybr2r1 = lsim( Mbr2, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
ybr2r2 = lsim( Mbr2, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
ybr2r3 = lsim( Mbr2, [ 0 0 1 ]' * ones( size( t ) ), t ) ;

% Figure 11.11

figure( 7 ) ; subplot( 2, 3, 1 ) ;
set( gcf, 'DefaultAxesColorOrder', [ 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1 ] ) ;
plot( t, yr1, '-', t, ybr2r1, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 2 ) ;
plot( t, yr2, '-', t, ybr2r2, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 3 ) ;
plot( t, yr3, '-', t, ybr2r3, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;


% 2. Balanced truncation.

ybt2r1 = lsim( Mbt2, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
ybt2r2 = lsim( Mbt2, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
ybt2r3 = lsim( Mbt2, [ 0 0 1 ]' * ones( size( t ) ), t ) ;

% Figure 11.12

figure( 8 ) ; subplot( 2, 3, 1 ) ;
set( gcf, 'DefaultAxesColorOrder', [ 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1 ] ) ;
plot( t, yr1, '-', t, ybt2r1, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 2 ) ;
plot( t, yr2, '-', t, ybt2r2, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 3 ) ;
plot( t, yr3, '-', t, ybt2r3, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;


% 4. Hankel norm approximation.

ybh2r1 = lsim( Mbh2, [ 1 0 0 ]' * ones( size( t ) ), t ) ;
ybh2r2 = lsim( Mbh2, [ 0 1 0 ]' * ones( size( t ) ), t ) ;
ybh2r3 = lsim( Mbh2, [ 0 0 1 ]' * ones( size( t ) ), t ) ;

%Figure 11.13

figure( 9 ) ; subplot( 2, 3, 1 ) ;
set( gcf, 'DefaultAxesColorOrder', [ 1 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1; 0 0 1 ] ) ;
plot( t, yr1, '-', t, ybh2r1, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 2 ) ;
plot( t, yr2, '-', t, ybh2r2, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;

subplot( 2, 3, 3 ) ;
plot( t, yr3, '-', t, ybh2r3, '--' ) ;
axis( SimStepAxis ) ;
xlabel( 'Time' ) ;