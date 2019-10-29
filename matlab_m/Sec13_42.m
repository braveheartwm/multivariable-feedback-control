%Section 12.4.2. Detailed LV-model
clear all;close all;
%Model (12.19)
% Distillation column model with 5 states and disturbances
% See Hovd and Skogestad, Auomatica, p. 989-996, 1992. 

A = [ -5.131e-3 0 0 0 0; 0 -7.366e-2 0 0 0; 0 0 -1.829e-1 0 0;
       0 0 0 -4.620e-1  9.895e-1; 0 0 0 -9.895e-1 -4.620e-1];

B = [-.629 .624; .055 -0.172; 0.030 -0.108; 
      -0.186 -0.139; -1.23 -0.056];


C= [-0.7223 -0.5170 0.3386 -0.1633e-1 0.1121;
    -0.8913 0.4728 0.9876 0.8425 0.2186];

D = [ 0 0; 0 0];


Bd = [-0.062 -0.067; 0.131 0.040; 0.022 -0.106; 
     -0.188 0.027; -0.045 0.014];

G = ss(A,B,C,D);
Gd = ss(A,Bd,C,D);
[Ad,Bd,Cd,Dd]=ssdata(Gd);

%-------------------------------------------------------------
%Analysis of model
%-------------------------------------------------------------

zero(G)

% Steady-state analysis
G0 = -C*inv(A)*B + D;
Gd0 = -Cd*inv(Ad)*Bd + Dd;
RGA0=vrga(G0);
PRGA0 = diag(diag(G0))/G0;
CLDG0 = PRGA0*Gd0;
GinvGd0 = inv(G0)*Gd0;


w = logspace(-3,1,41);
svdG = wsvd(G,w); Gdf=frd(Gd,w);
for i=1:length(w)
    Gf=freqresp(G,w(i));
    RGA_w(:,:,i) = Gf.*inv(Gf).';
end
RGA = frd(RGA_w,w);


%%%%%%%%%%%% MAKE FIGURE 13.16:
subplot(221)
loglog(w,svdG,w,1,':');
axis([0.001,10,.1,300]);drawnow
subplot(222)
uplot('liv,lm',RGA(1,1),RGA(1,2));
axis([0.001,10,.1,300]);




