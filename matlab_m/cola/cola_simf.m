% Program cola_simf.m
% Simple dynamic simulation of encrease in feedrate)
% You must first run: cola_init.m

disp('Simulating 1% increase in feed rate with all flows (L,V,D,B) constant....');
% FIRST GO INTO THE FILE cola.m and change F=1.00 to F=1.01. 
% SAVE THE FILE as cola_F1.m_and simulate:

[t,x]=ode15s('cola4_F1',[0 500],Xinit); 
t0 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting

disp('Finished. Plot MB, xB and yD as a function of time.....');
% Plot reboiler holdup (Apparent delay of about 1.26 min)
figure(2); plot(t0,M1); axis([0 3 0.5 0.52]); title('MB: reboiler holdup')

figure(3); plot(t0,xB); title('xB: reboiler composition')
figure(4); plot(t0,yD); title('yD: distillate composition')
