%Helicopter model (updated in 2005 with newer matlab commands)
%For more details with H-infinity design etc., see original file from first edition (1996): 
%http://www.nt.ntnu.no/users/skoge/book/1st_edition/matlab_m/Sec12_2.m

a01 = [          0                  0                  0   0.99857378005981;
                 0                  0   1.00000000000000  -0.00318221934140;
                 0                  0 -11.57049560546880  -2.54463768005371;
                 0                  0   0.43935656547546  -1.99818229675293;
                 0                  0  -2.04089546203613  -0.45899915695190;
-32.10360717773440                  0  -0.50335502624512   2.29785919189453;
  0.10216116905212  32.05783081054690  -2.34721755981445  -0.50361156463623;
 -1.91097259521484   1.71382904052734  -0.00400543212891  -0.05741119384766];

a02 = [0.05338427424431             0                  0                  0;
  0.05952465534210                  0                  0                  0;
 -0.06360262632370   0.10678052902222  -0.09491866827011   0.00710757449269;
                 0   0.01665188372135   0.01846204698086  -0.00118747074157;
 -0.73502779006958   0.01925575733185  -0.00459562242031   0.00212036073208;
                 0  -0.02121581137180  -0.02116791903973   0.01581159234047;
  0.83494758605957   0.02122657001019  -0.03787973523140   0.00035400385968;
                 0   0.01398963481188  -0.00090675335377  -0.29051351547241];

a0=[a01 a02];

b0=[              0                  0                  0                  0;
                  0                  0                  0                  0;
   0.12433505058289   0.08278584480286  -2.75247764587402  -0.01788876950741;
  -0.03635892271996   0.47509527206421   0.01429074257612                  0;
   0.30449151992798   0.01495801657438  -0.49651837348938  -0.20674192905426;
   0.28773546218872  -0.54450607299805  -0.01637935638428                  0;
  -0.01907348632812   0.01636743545532  -0.54453611373901   0.23484230041504;
  -4.82063293457031  -0.00038146972656                  0                 0];

c0 = [ 0        0         0         0         0    0.0595   0.05329  -0.9968;
     1.0        0         0         0         0         0         0        0;
       0      1.0         0         0         0         0         0        0;
       0        0         0  -0.05348       1.0         0         0        0;
       0        0       1.0         0         0         0         0        0;
       0        0         0       1.0         0         0         0       0];

d0 = zeros(6,4);

sys = ss(a0,b0,c0,d0);

%Remove fast states
tol = 2; % to eliminate poles at -11.49 and -2.3036
k = 6;
p = pole(sys);

sysd=canon(sys); % Diagonalize the system 
elim=(abs(p)>tol) & (real(p)<0); % fast stable states (abs(x) > tol)
syst=modred(sysd,elim,'t'); % then: Truncate fast states. 
sysr=modred(sysd,elim); % or: Residualize fast states.

%Balanced model reduction

n=size(sys.A,1); 
sysbt=balancmr(sys,k);     % kth order balanced truncation. 
%sysbt=modred(balreal(sys),k+1:n,'t'); % Alternate method for balanced truncation  
sysbr=modred(balreal(sys),k+1:n); % or: kth order balanced residualization.
sysbh=hankelmr(sys,k); % or: kth order optimal Hankel norm approx.

%Using coprime factors
nu=size(sys,2);  
sysct=ncfmr(sys,k); % balanced truncation of coprime factors. 
[sysc,cinfo]=ncfmr(sys,n); % or: obtain coprime factors 
syscr=modred(cinfo.GL,k+1:n); % and residualize. 
syscrm=minreal(inv(syscr(:,nu+1:end))*syscr(:,1:nu)); % and obtain kth order model. 
sysch=hankelmr(cinfo.GL,k); % or: optimal Hankel norm approximation. 
syschm=minreal(inv(sysch(:,nu+1:end))*sysch(:,1:nu)); % obtain reduced order model. 
 
