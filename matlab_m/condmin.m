function gmin = condmin(Gf)
% Minimize cond.no. with both input and output scaling
% WARNING: May be inaccurate because of inaccurate mu-calculation
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: condmin.m,v 1.1 2004/01/22 08:19:04 vidaral Exp $

[n,m]=size(Gf);
zero = 0*eye(n);
H = [zero inv(Gf); Gf zero];
blk0=[1 1]; blk1=blk0;
for k=2:n
  blk1 = [blk1; blk0];
end
blk2=blk1; blk=[blk1; blk2];
mu1=mussv(H,blk,'c'); gmin1=mu1(:,1);
gmin=gmin1*gmin1;
