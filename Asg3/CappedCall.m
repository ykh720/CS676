function V = CappedCall(S0, r, sigma, pu,mu_u,mu_d, lambda, K,T,C,M,N)
% we add M number of simulations and N number of timesteps as input
% parameters

dt = T/ N;

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

% obtain MC capped call option values
V = mean(exp(-r*T) * min(max(S-K,0),C));
end

