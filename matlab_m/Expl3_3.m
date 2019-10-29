%Example 3.3
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl3_3.m,v 1.1 2004/01/22 13:55:32 standber Exp $

d10=1;
points = 100;
d20 = linspace(-5,5,points);
d=[d10*ones(1,points);d20];
G=[ 5 4; 3 2];
y=G*d;
d_norm=sqrt(d(1,:).*d(1,:)+d(2,:).*d(2,:));
y_norm=sqrt(y(1,:).*y(1,:)+y(2,:).*y(2,:));
out=y_norm./d_norm;
figure(1);
clf;
plot(d20,out)
xlabel('d20 / d10');
ylabel('|| y ||_2  / || d ||_2');
text(-4.5,7.34,'MAXIMUM  SINGULAR VALUE = 7.34')
text(-0.7,0.37,'MINIMUM SINGULAR VALUE = 0.27')

