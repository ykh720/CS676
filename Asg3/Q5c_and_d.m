%% Q5c and d
clearvars
close all
% Parameters 
rng('default')

C = Inf;
sigma = 0.15;
r = 0.05; % we use r, but should we consider it as compensated?
T = 1;
Klist = linspace(70,120,20);
S0 = 95;
mu_u = 3.04;
mu_d = 3.08;
pu = 0.34;
lambda = 0.1;
N = 800;
M = 25000;

V = zeros(length(Klist), 1);

for i = 1:length(V)
    K = Klist(i);
    V(i,1) = CappedCall(S0, r, sigma, pu,mu_u,mu_d, lambda, K,T,C,M,N);
end

V

IV = blsimpv(S0, Klist', r, T, V); 
% note that default class for blsimpv is call

g = figure(1);
plot(Klist, IV')
title('Implied Volatility against Strike')
xlabel('Strike')
ylabel('Implied Volatility')
saveas(g, 'q5c','epsc')


%% Q5d
% We see a volatility skew, that is implied volatility decreases as the
% strike increases. We know that for equity options, the implied volatility
% against strike is often downward sloping. So this means that the assumed
% jump model does model such observed real life phenomenon. This may
% suggest that it is a good model to model the market underlying movement.
% Note that if we just use the geometric brownian motion with no jumps, we
% are just obtaining approximation of Black-Scholes price and the implied
% volatility plot will be close to a constant.
% 
% There are also some oscillations when strike is large but I think it is
% due to the numerical approximation error (sampling error and time
% discretization error).
%
% Note also that the calculated implied volatilities are all larger than
% 0.25, which is again large than simga = 0.15. This means that our jump
% model "adds" more volatility (in BS sense) via allowing the underlying
% price to have jumps. This is also expected.