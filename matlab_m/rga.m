function lambda=rga(g)
%   Lambda=RGA(G) 
%
%   Calculate the relative gain array (RGA) of a matrix G.
%
%   By Yi Cao 02/03/95
%

lambda=g.*(pinv(g)).';
