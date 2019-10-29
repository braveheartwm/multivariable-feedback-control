
% This is file cola_commands.m
% It contains are some useful commands for MIMO controllability analysis:
% Poles/zeros and their directions, RGA, PRGA, CLDG, singular values
% Uses the Mu-toolbox.

% Start from G and Gd: 

% start model:
% Here: Distillation column model with 5 states and disturbances
% Scaled LV-model (Example 10.10 in book)
% The two inputs are (LV); 
% the two disturbances are (F,zF)
% The two outputs are (yd,xb)
A = [ -5.131e-3 0 0 0 0; 0 -7.366e-2 0 0 0; 0 0 -1.829e-1 0 0;
       0 0 0 -4.620e-1  9.895e-1; 0 0 0 -9.895e-1 -4.620e-1]
B = [-.629 .624; .055 -0.172; 0.030 -0.108; 
      -0.186 -0.139; -1.23 -0.056]
C= [-0.7223 -0.5170 0.3386 -0.1633e-1 0.1121;
    -0.8913 0.4728 0.9876 0.8425 0.2186]
D = [ 0 0; 0 0]
Bd = [-0.062 -0.067; 0.131 0.040; 0.022 -0.106; 
     -0.188 0.027; -0.045 0.014]
G = pck(A,B,C,D); Gd = pck(A,Bd,C,D);
% end model

% We now have G and Gd .....  Start analysis.
format short e

spoles(G)                        % Compute poles   (Distillation: "slow" is at 0.0052=1/194 min)
[A,B,C,D]=unpck(G); poles=eig(A) %   should be the same ...
[T,Po] = eig(A);  YP = C*T       %   pole output vectors (to obtain directions normalize columns to 1)
[Q,Pi] = eig(A'); Q=conj(Q);
UP = B'*Q                        %   pole input vectors
Shouldbezero=Po-Pi;               % if nonzero: WARNING - ordering of eigenvalue differs in YP and UP
if norm(Shouldbezero)>1.e-7 disp('WARNING: Eigenvalue order differs for inputs and output pole vectors'),
    disp ('Use opde and ipde instead:  [Po, Yp, Xpo, Spo] = opde(G), [Pi, Up, Xpi, Spi] = ipde(G)'); end

Szeros=szeros(G)                 % zeros (the ones "far away" are not important)
Tzero=tzero(A,B,C,D)            %   should be the same ... (If they differ I would trust tzero)
[Zy,YZ,Xy] = ozde(G); Zy, YZ     %   zero output directions
[Zu,UZ,Xu] = izde(G); Zu, UZ     %   zero input directions       
           
G0 = frsp(G,0)                   % steady-state plant gains
G00= -C*inv(A)*B + D             %   should be the same ...
RGA0 = G00.*pinv(G00.')          % steady-state RGA
vrga(G0)                         %   should be the same ... (using subroutine vrga)
PRGA0 = mmult(vdiag(vdiag(G0)),minv(G0)) % Performance RGA
Gd0= frsp(Gd,0)                  % steady-state disturbance gains
CLDG0 = mmult(PRGA0,Gd0)         % Closed-loop disturbance gain
GinvGd0 = mmult(minv(G0),Gd0)    % inputs for perfect disturbance rejection

w = logspace(-3,1,41);           % Now do the same as a function of frequency
Gf = frsp(G,w); Gdf=frsp(Gd,w);  % Frequency repsonse of G and Gd
vplot('liv,lm',vsvd(Gf),1,':');  % Plot singular values of G
vplot('liv,lm',Gdf,1,':');       % Plot elements in Gd
vplot('liv,lm',vrga(Gf),1,':');  % RGA-elements
gdiag=vdiag(vdiag(Gf)); prga = mmult(gdiag,minv(Gf));
vplot('liv,lm',prga,1,':');      % PRGA-elements
cldg=mmult(prga,Gdf);
vplot('liv,lm',cldg,1,':');      % CLDG-elements

% Some closed-loop analysis ..........

% Two simple PI controllers
k1 = nd2sys([3.76 1], [3.76 1.e-4],  0.261);
k2 = nd2sys([3.31 1], [3.31 1.e-4], -0.375);
K = daug(k1,k2); GK = mmult(G,K);
S = minv(madd(eye(2),GK));              % sensitivity function
spoles(S)                               % closed-loop poles (stable)
Sf=frsp(S,w); vplot('liv,lm',Sf,1,':'); % sensitivity function elements
vplot('liv,lm',vsvd(Sf),1,':');         % singular values of sensitivity 
SGd = mmult(S,Gd); SGdf=frsp(SGd,w);    
vplot('liv,lm',vsvd(SGdf),1,':');  % closed-loop effect of disturbances 
                                   % (which is less than 1 as desired)
y = trsp(SGd,[1 0]',100,0.2); vplot(y); % Time response to disturbance in F 
y = trsp(SGd,[0 1]',100,0.2); vplot(y); % Time response to disturbance in zF 


