clear all; close all

s = tf('s');
G =  (0.01/(s+1.72e-4)/(4.32*s + 1))*[-35.54*(s+0.0572), 1.913; -30.22e5*s, -9.188e5*s-638.56];
omega = logspace(-5,2,61);

for i = 1:length(omega)
    Gf = freqresp(G,omega(i));
    RGA_w(:,:,i) = Gf.*inv(Gf).';
    RGAnumber_diag(i) = sum(sum(abs(RGA_w(:,:,i) - eye(2)))); 
    RGAnumber_offdiag(i) = sum(sum(abs(RGA_w(:,:,i) - [0 1;1 0]))); 
end
RGA = frd(RGA_w,omega);

%Plotting magnitudes of RGA elements
% For 2 x 2 systems, RGA(2,1) = RGA(1,2) and RGA(2,2) = RGA(1,1)
plot(abs(RGA(1,1)))
hold
plot(abs(RGA(1,2)),'r')
axis([1e-5 100 0 1])
text(0.08,0.32,'|L_{11}| = |L_{22}|'),text(0.27,0.52,'|L_{12}| = |L_{21}|')
ylabel('|L_{ij}|'), xlabel('\omega')

%Plotting RGA number
figure, 
plot(omega, RGAnumber_diag)

hold
plot(omega, RGAnumber_offdiag,'r:')
text(0.07,1.3,'Off-diagonal pairing'),text(0.32,2,'Diagonal pairing')
ylabel('||L - I||_{sum}'), xlabel('\omega')