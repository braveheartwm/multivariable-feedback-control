



function SVD = wsvd(G,w)
%WSVD SVD = WSVD(G,W) returns frequency varying maximum and minimum SVDs 
%
% $Id: wsvd.m,v 1.1 2005/11/01 09:00 hori Exp $

if (nargin == 0) | (nargin > 2),
   disp('usage: wsvd(G,w)')
   return
end

for i=1:length(w)
    Gf=freqresp(G,w(i));
    SVD(i,1)=max(svd(Gf));
    SVD(i,2)=min(svd(Gf));
end

