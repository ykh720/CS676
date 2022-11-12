%% Q5b
% Note that we haven't used CappedCall function this time. Instead we
% basically implement the function in this file and compute the fair values
% of the capped call. The way we used here will be slightly faster since we
% only simulate once and the resulting paths will be used across all values
% of cap C. Of course, we could do it by calling CappedCall for each C. 

% Parameters 
close all

rng('default')

sigma = 0.15;
r = 0.05; % we use r, but should we consider it as compensated? I think so.
T = 1;
K = 95;
S0 = 95;
mu_u = 3.04;
mu_d = 3.08;
pu = 0.34;
lambda = 0.1;
% dt = 1/1000;    
% use 1/1000 or N = 800?

N = 800;
dt = T/ N;
M = 25000;
Clist = 20:10:100;

% compensated drift E[J-1]
kappa = pu* mu_u / (mu_u - 1) + (1-pu) *mu_d /(mu_d +1) -1;
% compensated drift for X = log(S), risk neutral
drift = (r - sigma^2/ 2 - lambda* kappa);

% X = log(S)
X_old = log(S0) * ones(M,1 );
X_new = zeros(M,1);

jump_check = zeros(M,1);
jump_size = zeros(M,1);
jump_mask = zeros(M,1);
jump_up = zeros(M,1);
jump_down = zeros(M,1);

% apply Euler timestepping on X = log(S)
for i = 1:N %timestep loop
    % uniform distribution 
    jump_check = rand(M,1);
    jump_mask = jump_check <= lambda *dt; % index for existence of jumping
    % resample, now for determining up or down jump
    jump_check = rand(M,1);
    jump_up = jump_check <= pu; % storing indices for up jump
    jump_down = ~jump_up;
    jump_size(jump_up) = exprnd(1/mu_u, sum(jump_up),1);
    jump_size(jump_down) = -exprnd(1/mu_d, sum(jump_down),1);

    jump_size = jump_size .* jump_mask;

    X_new = X_old + drift* dt + sigma *sqrt(dt) * randn(M,1) + ...
        jump_size;
    X_old = X_new;
end

S = exp(X_new); % asset price values for each path

% do we have to generate different samples for different C? I think we
% don't. But of course we can easily do it.

% obtain MC capped call option values
V = zeros(length(Clist),1);
for i = 1:length(Clist)
    C = Clist(i);
    V(i,1) = mean(exp(-r*T) * min(max(S-K,0),C));
end

V_mat = V;
V = array2table(V, 'RowNames', "C = " + string(Clist))

g = figure(1);
plot(Clist, V_mat)
title('Capped call prices against cap C')
xlabel('C')
ylabel('Capped call prices')
saveas(g, 'q5b','epsc')

% V_fun = CappedCall(S0, r, sigma, pu,mu_u,mu_d, lambda, K,T,Clist(1),M,N)
% 
% for i = 1 : length(Clist)
%     V_fun_list(i) = CappedCall(S0, r, sigma, pu,mu_u,mu_d, lambda, K,T,Clist(i),M,N);
% end
% 
% V_fun_list 

%% Discussion
% How does the computed option value depend on the cap C ? Explain why your
% observation is reasonable.
%
% We see that the computed option value increases as the cap C increases.
% This is reasonable since the cap C is limiting the highest possible final
% payoff we can get from the capped call. If S_T is very large, let say S_T
% - K = 10000, then a capped call with C = 20 will only have final payoff
% 20 while a capped call with C = 100 will have final payoff 100 instead.
% In other words, a higher cap will allow us to earn more on some extreme
% in the money asset path so that the resulting option value is higher.
%
% Also note that the speed of increasing for option value as C increases is
% actually decreasing. This is reasonable since as the cap C increases, the
% number of asset price path that can exceed such cap becomes rarer (in
% fact very rare, the number of asset price path exceeding certain C,
% denoted by g(C), perhaps of order of an inverse of high order polynomial
% in C, or even inverse of an exponential function in C. I think this can
% be investigated further by examining brownian motion properties).
