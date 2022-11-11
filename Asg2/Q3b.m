%% Parameters initialization 

K = 102;
B = 100;
T = 0.5; % expiry of 6 months
x = 15;
sigma = 0.2;
r = 0.03;
S0 = 105;

% Discussion for time discretization error
% 
% The error in the computed  $\Tilde{V}$ actually depend on the time
% discretization. Note that barrier option is a path-dependent option,
% meaning that its payoff not only depend on the final value of the stock
% but also the path taken to reach this final value. The way that we
% approximate the payoff is based on a time discretization. That is, if our
% simulation suggests that at each time step $t_n$, we have $S(t_n) > B$
% and also we have $S(T) < K$, then we receive the cash amount x. However,
% between successive time step $t_n, t_{n+1}$, although we have $S(t_n),
% S(t_{n+1}) >B$, there is no guarantee that $S(s) > B$, for all $s$
% strictly between $t_n$ and $t_{n+1}$. If this is the case, then by the
% definition of down barrier option, we cannot receive cash amount $x$ in
% the final time. This produces an error where by our simulation we
% approximate the payoff by $x$ whereas the true payoff is actually 0. We
% see that the computed  $\Tilde{V}$ actually depend on the time
% discretization, in the sense that the time discretization fails to
% capture the possible fluctuation of price during the time between
% successive time discretization steps.

% Code up the MC algorithm

dtlist = [1/200, 1/400, 1/800, 1/1600, 1/3200, 1/6400];
Mlist = 1000:3000:100000; 

rng('default') % for reproducibility 

output = zeros(length(dtlist), length(Mlist));
for i = 1: length(dtlist)
    for j = 1:length(Mlist)
        output(i,j) = MC_barrier(sigma, r, T, K, S0, x, B, dtlist(i), Mlist(j));
    end
end

% obtain the exact price by part a
exact = do_cashornth_put(sigma, r, T, K, S0, x, B);

% for each of Delta t, plot the computed option values using MC for
% M=1000:3000:100000
for i = 1:length(dtlist)
    h(i) = figure(i);
    plot(Mlist, output(i, :))
    xlabel('number of sample paths M')
    ylabel('option value')
    title(strcat('timestep: ',string(dtlist(i))))
    hold on
    plot(Mlist, exact* ones(length(Mlist),1) ,'r-*')
    legend('MC price','exact')
    hold off
end

% We observe from the plots that as M goes up, that is more path
% realizations are used, the MC option value is closer to the exact price.
% However, such improvements are not very significant in the case of large
% timesteps like 1/200, 1/400, 1/800. It is more significant for the
% small timestep like 1/6400. It is because in such case the time
% discretization error (presumably O(Delta t)) is small and so the dominant
% error is sampling error O(1/sqrt(M)). In order to match these two errors
% and achieve optimal error, we need M to be O(1/Delta t ^ 2). In the case of
% timestep 1/6400, since 1/(1/6400)^2 = 40960000, way larger than 100000,
% increasing M will give us more significant improvement. Another
% observation is that when M goes up, the MC values seem to stabilize. This
% can intuitively explained by that increase in M will reduce the effect of
% random path realization.

% Here are some other plots that I think might be of interest

%% Explanation for figure g(1)
% Here we fix the number of sample paths to be 100000 and plot the option
% value against different time discretization. We see the general pattern
% that the finest the time discretization is, the more accurate the MC
% price is. This is becasuse we have time discretization error and using a
% finer time discretization scheme will reduce such error.

g(1) = figure;
plot(dtlist, output(:,end))
xlabel('Delta t')
ylabel('option value')
title('100000 sample paths for each time discretization')
hold on 
plot(dtlist, exact* ones(length(dtlist),1) ,'r-*')
legend('MC price','exact')
hold off
% Consider the smallest (or finest) time discretization scheme and plot the
% MC price for different number of sample paths

%% Explanation of figure g(2)
% Here we plot a 3D plot of the absolute value of the difference between MC
% price and exact price. We observe that as time discretization becomes
% finer and finer, the accuracy of MC price increase significantly. This is
% because time discretization is reduced. 
% Moreover, we see that as M increases, the MC price becomes stabler. This
% is because the reduction of the effect of randomness. We also see slight
% increase in MC price accuracy as M increases, due to reduction in
% sampling error.
% 
% The value of the sampling error and time discretization error can be
% somewhat estimated by observing the following 
% 1./ sqrt(Mlist)
% dtlist

g(2) = figure;
[X,Y] = meshgrid(Mlist, dtlist);
mesh(X,Y,abs(output- exact))
ylabel('Delta t')
xlabel('M')
title('3D plot of absolute error')

% Save the figures
for i=1:length(h)
    saveas(h(i),strcat('fig_needed',string(i)),'epsc')
end
for i =1:length(g)
    saveas(g(i),strcat('fig_extra3',string(i)),'epsc')
end

function value = MC_barrier(sigma, r, T, K, S0, x, B, dt, M)
% MC_barrier returns the down-and-out cash-or-nothing put option
%   value by Monte Carlo simulation
% 
% Input: 
% sigma: volatility of the underlying 
% r: interest rate
% T: time of expiry
% K: strike price
% S0: initial asset value
% x: cash payout
% B: down barrier 
% dt: timestep
% M: number of path realization 

N = T/ dt; % assume it is integer
V = zeros(M,1); % placeholder for the option value of each path
idx = true(M,1); % storing the indices of path that has positive final payout 
S = S0 * ones(M,1); % initial asset value 

for n = 0:N-1
    % obtain phi_n for each path
    sample = randn(M,1);
    % calculate S(t_{n+1}), see (8) in assignment
    S = S.*exp((r - 1/2 * sigma^2) * dt + sigma.* sample .* sqrt(dt));
    % record the indices of path of S that is below the barrier at this timestep 
    check = S <= B;
    idx(check) = false;
end

% final check: only realization with final asset value less than K can have
% cash payout
endcheck = S >= K;
idx(endcheck) = false;

% compute the MC price by (9) in assignment 
V(idx) = x;
value = exp(-r*T) *mean(V);
end 

function V = do_cashornth_put(sigma, r, T, K, S0, x, B)
% do_cashornth_put returns the down-and-out cash-or-nothing put option
%   value
%
% allow vector input
%
% Input: 
% sigma: volatility of the underlying 
% r: interest rate
% T: time of expiry
% K: strike price
% S0: initial asset value
% x: cash payout
% B: down barrier 

z_1 = log(S0./K)./ (sigma .* sqrt(T)) + sigma.*sqrt(T)/2;
z_2 = log(S0./B)./ (sigma .* sqrt(T)) + sigma.*sqrt(T)/2;
y_1 = log(B.^2./(S0.*K))./ (sigma .* sqrt(T)) + sigma.*sqrt(T)/2;
y_2 = log(B./S0)./ (sigma .* sqrt(T)) + sigma.*sqrt(T)/2;

inter = sigma.*sqrt(T);

V = x.*exp(-r.*T) .*(normcdf(-z_1 + inter) - normcdf(-z_2 + inter) ...
    + (S0./B) .* normcdf(y_1 - inter) - (S0./B).*normcdf(y_2 - inter));
end

