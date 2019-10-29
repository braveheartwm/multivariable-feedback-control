% Section 13.3.2          Aero-engine control, Part I
%
% This file is used to produce the plots regarding analysis of output
% selection problem for the aero engine.
% 
% Dependency: data files aero1.mat, aero2.mat
%
% Kjetil Havre, 19/12-1995.
% Modified Federico Zenith, 2004
% Modified Vidar Alstad, 2004
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
%

% $Id: Sec13_32.m,v 1.7 2004/04/20 13:35:13 vidaral Exp $

clear all;close all;

% Load the engine model
load aero1
w = logspace(-3,3,301);			

% ********** FOLLOWING LINES MAY BE USELESS (NO Si or So PRESENT)
% Scale the plant G_eng (Si and So are the input & output scaling matrices)
% G = So * G_eng * Si ;
G=Gall;
Gjw=freqresp(G,w);
G0=freqresp(G,0);
[U0,S0,V0]=svd(G0);

 
Goh=round(1000*G0(1:6,1:3))/1000;
rgaGoh=Goh.*pinv(Goh)';

% Displaying steady state and RGA

disp('Gall(0)=');
fprintf(1, '     %6.3f  %6.3f  %6.3f\n', Goh' )

disp('RGA(Gall)=');
fprintf( 1, '     %6.3f  %6.3f  %6.3f\n', rgaGoh' )

disp('RGA_sum=');
fprintf( 1, '     %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n', sum( rgaGoh' ) )


% Displaying singular value decomposition

U0h = round( 1000 * U0( 1:6, 1:6 ) ) / 1000;
disp( 'U0=' );
fprintf( 1, '     %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n', U0h' )

S0h = round( 1000 * S0( 1:6, 1:3 ) ) / 1000;
disp( 'S0=' ) ;
fprintf( 1, '     %6.3f  %6.3f  %6.3f\n', S0h' )

V0h = round( 1000 * V0( 1:3, 1:3 ) ) / 1000;
disp( 'V0=' );
fprintf( 1, '     %6.3f  %6.3f  %6.3f\n', V0h' )

% Zeros

Z  =tzero(G);
Z1 =tzero(G(S1,:)); 
Z2 =tzero(G(S2,:));
Z3 =tzero(G(S3,:));
Z4 =tzero(G(S4,:));
Z5 =tzero(G(S5,:));
Z6 =tzero(G(S6,:));


% We have RHP zeros < 100 rad/s for the following sets:
%    In set 3, z = 30.853 ==> not use both NL and LPEMN measurements
%    In set 6, z = 27.664 ==> not use both OPR2 and LPEMN measurements
%
%    NL, LP-compressor spool speed is a function of LPPR the same is 
%    LPEMN ++> interacting effects which results in RHP zero
%
%    The OPR2 is a function of both LPPR and HPPR. LPEMN is a function of
%    LPPR ++> interacting effects, which results in RHP zero.
%
%    Conclusion: do not need to consider set S3 and S6
% clear Z3 Z6 s3 s6


disp( ' ' ) ;
disp( '----------------------------------------' );
disp( 'Frequency-dependent relative gain array.' );
disp( '----------------------------------------' );

disp( 'Frequency-dependent RGA for full process.' );

RGA10 =G0(S1([2,1,3]),:).*pinv(G0(S1([2,1,3]),:))'; 
RGA20 =G0(S2([2,1,3]),:).*pinv(G0(S2([2,1,3]),:))';
RGA30 =G0(S3([2,1,3]),:).*pinv(G0(S3([2,1,3]),:))';
RGA40 =G0(S4,:).*pinv(G0(S4,:))';
RGA50 =G0(S5,:).*pinv(G0(S5,:))';
RGA60 =G0(S6,:).*pinv(G0(S6,:))';

for i = 1:size(Gjw,3)
    RGA(:,:,i)=Gjw(:,:,i).*pinv(Gjw(:,:,i)');
   RGA1(:,:,i)=Gjw(S1([2,1,3]),:,i).*pinv(Gjw(S1([2,1,3]),:,i).');
   RGA2(:,:,i)=Gjw(S2([2,1,3]),:,i).*pinv(Gjw(S2([2,1,3]),:,i).');
   RGA3(:,:,i)=Gjw(S3([2,1,3]),:,i).*pinv(Gjw(S3([2,1,3]),:,i).');
   RGA4(:,:,i)=Gjw(S4,:,i).*pinv(Gjw(S4,:,i).');
   RGA5(:,:,i)=Gjw(S5,:,i).*pinv(Gjw(S5,:,i).');
   RGA6(:,:,i)=Gjw(S6,:,i).*pinv(Gjw(S6,:,i).');
end

% Calculation of frequency-dependent RGA numbers

for i=1:size(Gjw,3)
    rn1(i)=sum(sum(abs(RGA1(:,:,i)-eye(3))));
    rn2(i)=sum(sum(abs(RGA2(:,:,i)-eye(3))));
    rn3(i)=sum(sum(abs(RGA3(:,:,i)-eye(3))));
    rn4(i)=sum(sum(abs(RGA4(:,:,i)-eye(3))));
    rn5(i)=sum(sum(abs(RGA5(:,:,i)-eye(3))));
    rn6(i)=sum(sum(abs(RGA6(:,:,i)-eye(3))));
end

% Figure 13.10

figure(1);
semilogx(w,rn1,w,rn2,w,rn3,w,rn4,'--',w,rn5,'-.',w,rn6);
axis([1e-3 1e2 0 7]);
xlabel('Frequency [ rad/s ]');
ylabel('RGA number');
legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'TL');

% Smallest singular values:

Smsv1=min(svd(G0(S1,:)));
Smsv2=min(svd(G0(S2,:)));
Smsv3=min(svd(G0(S3,:)));
Smsv4=min(svd(G0(S4,:)));
Smsv5=min(svd(G0(S5,:)));
Smsv6=min(svd(G0(S6,:)));

fprintf( '\nThe smallest steady-state singular value for set 1 is %0.3f;', Smsv1 );
fprintf( '\nThe smallest steady-state singular value for set 2 is %0.3f;', Smsv2 );
fprintf( '\nThe smallest steady-state singular value for set 3 is %0.3f;', Smsv3 );
fprintf( '\nThe smallest steady-state singular value for set 4 is %0.3f;', Smsv4 );
fprintf( '\nThe smallest steady-state singular value for set 5 is %0.3f;', Smsv5 );
fprintf( '\nThe smallest steady-state singular value for set 6 is %0.3f.\n\n', Smsv6);


load aero2

Gs=So*G_eng*Si;

P =   [ 0 0 0 1 0 0;
        1 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 1 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0];

S4=[ 1 5 6 ];
S5=[ 2 5 6 ];

Gall=P*Gs;
G4=Gall(S4,:);
G5=Gall(S5,:);

G1i=G4(1,:);      % Output y1 = NL
G2i=G5(1,:);      % Output y2 = OPR1
[G2ib,H4]=balreal(G1i);
[G2ib,H5]=balreal(G2i);

k=1:1:15;


% Figure 13.11

figure(2);
plot(k,H4,':',k,H4,'+',k,H5,':',k,H5,'x');
axis([0 16 0 0.4]);
text(10, 0.3, 'Set 4: +');
text(10, 0.2, 'Set 5: x');
xlabel( 'n-state number');
