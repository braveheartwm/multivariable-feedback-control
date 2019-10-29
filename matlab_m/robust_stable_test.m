clear all;
close all;
clc;

% 频率范围
omega = logspace(-3,2,101);

% 控制对象
tau = 75;
G = tf([0 1],[tau 1])*[-87.8 1.4;-108.2 -1.4];

% 控制器
K = tf([tau 1],[1 0])*[-0.0015 0;0 -0.075];

Wi=tf([1 0.2],[0.5 1]);

Wif = frd(Wi,omega);

LI = G*K;
SI = inv(eye(2)+LI);
TI = SI*LI;
Tjw = sigma(TI,omega);% 这里求取对应频率下的奇异值。

blk = [1 1;1 1];%该参数为mussv函数中的structure

[bnds,muinfo] = mussv(frd(TI,omega),blk);
bodemag(1/Wif,'r--',bnds(1,1),'-',omega)
hold on;
plot(omega,Tjw(1,:));
text( 1, 12, 'sv' ); text(0.002, 8, '1/WI' ); text( 1, 0.3, 'mu' );
hold off
grid;
axis auto