%Example 6.12
%Triangular plant
% Uses the functions rga, condmin, condmini, wcsens
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl6_9.m,v 1.3 2004/04/14 08:13:34 vidaral Exp $

clear all; close all;
G = [1 100; 0 1];
COND=cond(G)  %10^4
RGA=rga(G)   %identity
CONDMIN=condmin(G) %1
CONDMINI=condmini(G) %200

%Inverse based controller
c=1;
dynk=tf(c*1,[1 1.e-6]); 
Kinv=dynk*inv(G);

%With full block uncertainty
Gunc=G*(eye(2)+0.2*ultidyn('DeltaI',[2 2]));
Lunc=Gunc*Kinv;
wcs=wcsens(Lunc,eye(2),'So');
wcs_fb=wcs.So.MaximumGain.UpperBound

%With diagonal block uncertainty
Gunc=G*(eye(2)+[0.2*ultidyn('DeltaI1',[1 1]) 0; 0 0.2*ultidyn('DeltaI2',[1 1])]) ;
Lunc=Gunc*Kinv;
wcs=wcsens(Lunc,eye(2),'So');
wcs_db=wcs.So.MaximumGain.UpperBound

%Results:
 clc
 disp(sprintf('Condition number of G: %7.1f',COND)); %10^4
 disp(sprintf('RGA(G)=[%5.1f%5.1f]\n       [%5.1f%5.1f]',RGA')); %identity
 disp(sprintf('Minimized condition number: %5.1f',CONDMIN)); %1
 disp(sprintf('Input-minimized condition number: %6.1f',CONDMINI)); %200
 disp(sprintf('With full block uncertainty %8.2f',wcs_fb))
 disp(sprintf('With diagonal block uncertainty %8.2f',wcs_db))
