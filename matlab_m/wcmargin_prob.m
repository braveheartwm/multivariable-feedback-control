G = 3*tf([-2 1],conv([5 1],[10 1])); % Plant
Wi = tf([10 0.33],[10/5.25 1]);  % Uncertainty Weight
Delta = ultidyn('Delta', [1 1]); % Uncertainty
Gp = G * (1 + Wi*Delta); % Perturbed plant
K = 1.13*tf([12.7 1],[12.7 0]); % Nominal Controller

nom = wcmargin(Gp*K*gain); % no warnings here

% Want to find, by what factor, the gain can be reduced to have robust
% stability

gain = 0.1; % Since only stable systems are handled, make system stable and calculate

% Iterate on upper bound on Worst-Case gain margin
for i = 1:5
    x = wcmargin(Gp*K*gain);
    gain_iter = x.GainMargin(2);
    gain = gain*gain_iter
end

% final gain margin = 0.2852 (Correct solution 0.2743, trial and error)

norm(minreal(Wi*G*K*gain*inv(1+G*K*gain)),inf,1e-6) % answer 1.0475
    