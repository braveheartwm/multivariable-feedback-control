w = 0:0.01:2*pi;

d = [cos(w)' sin(w)']

G = [5 4; 3 2];

y = G * d';

figure('position',[150,200,900,380])
subplot(1,2,1)
plot(d(:,1), d(:,2))

subplot(1,2,2)
plot(y(1,:), y(2,:))

gain = sqrt(sum(y.^2));

[gain1,index1] = max(gain);
[gain2,index2] = min(gain);

hold
axis([-8,8,-4,4])
line([0 y(1,index1)],[0 y(2,index1)],'linestyle',':','linewidth',1)
line([0 y(1,index2)],[0 y(2,index2)],'linestyle',':','linewidth',1)
text(y(1,index1)+0.2, y(2,index1) + 0.1, 's1')
text(y(1,index2)-0.7, y(2,index2) + 0.25, 'sn')

subplot(1,2,1)
hold
line([0 d(index1,1)],[0 d(index1,2)],'linestyle',':','linewidth',1)
line([0 d(index2,1)],[0 d(index2,2)],'linestyle',':','linewidth',1)
text(0.2, 0.35, 'v1'), text(0.2, -0.5, 'vn')