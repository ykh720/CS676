%% Q3e

%% Data preparation 
clearvars
close all

rng('default') % reproducibility 

% parameter initialization

M = 50000; % as suggested on Piazza, suggested 2000
sigma = 0.2;
r = 0.03;
mu = 0.15;
T = 1;
S0 = 10;
K = 0.95 * S0;
N = 250;
dt = T/N;

[V0,S, delta_bin] = binomialDeltaStraddle(S0,r,sigma, T,N,K);

% no hedging
% Initial portfolio
S_MC = S0 * ones(M,1); 
B = (V0-delta_bin(1,1)*S0)*ones(M,1);
delta_hedge_old = delta_bin(1,1)*ones(M,1);
delta_hedge_new = delta_hedge_old;

% at t_N
sample = randn(M,1);
% calculate S(t_{n+1}), see (8) in assignment 2
S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * T + sigma.* sample .* sqrt(T));
W = max(S_MC -K, 0) + max(K- S_MC, 0);
hedgingerror = -W+ delta_hedge_new .* S_MC + B.*exp(r*T);

PL_nohedging = exp(-r*T) * hedgingerror ./ V0;

% daily hedge
% Initial portfolio
S_MC = S0 * ones(M,1); 
B = (V0-delta_bin(1,1)*S0)*ones(M,1);
delta_hedge_old = delta_bin(1,1)*ones(M,1);
delta_hedge_new = delta_hedge_old;
for n = 1:N-1
    % at time t_n
    % obtain phi_n for each path
    sample = randn(M,1);
    % calculate S(t_{n+1}), see (8) in assignment 2
    S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt + sigma.* sample .* sqrt(dt));
    delta_hedge_new = interpDelta(delta_bin(1:n+1, n+1), S(1:n+1, n+1), S_MC);
    B = B.*exp(r * dt) + (delta_hedge_old - delta_hedge_new).* S_MC;
    delta_hedge_old = delta_hedge_new;
end

% at t_N
sample = randn(M,1);
% calculate S(t_{n+1}), see (8) in assignment 2
S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt + sigma.* sample .* sqrt(dt));
W = max(S_MC -K, 0) + max(K- S_MC, 0);
hedgingerror = -W+ delta_hedge_new .* S_MC + B.*exp(r*dt);

PL_daily = exp(-r*T) * hedgingerror ./ V0;

% weekly hedge
% Initial portfolio
S_MC = S0 * ones(M,1); 
B = (V0-delta_bin(1,1)*S0)*ones(M,1);
delta_hedge_old = delta_bin(1,1)*ones(M,1);
delta_hedge_new = delta_hedge_old;
dt_week = 5*dt;
for n = 5:5:N-5
    % at time t_n
    % obtain phi_n for each path
    sample = randn(M,1);
    % calculate S(t_{n+1}), see (8) in assignment 2
    S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt_week + sigma.* sample .* sqrt(dt_week));
    delta_hedge_new = interpDelta(delta_bin(1:n+1, n+1), S(1:n+1, n+1), S_MC);
    B = B.*exp(r * dt_week) + (delta_hedge_old - delta_hedge_new).* S_MC;
    delta_hedge_old = delta_hedge_new;
end
% after the for loop, we get t_{N-5} hedging position and bond value

% at t_N
sample = randn(M,1);
% calculate S(t_{n+1}), see (8) in assignment 2
S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt_week + sigma.* sample .* sqrt(dt_week));
W = max(S_MC -K, 0) + max(K- S_MC, 0);
hedgingerror = -W+ delta_hedge_new .* S_MC + B.*exp(r*dt_week);

PL_weekly = exp(-r*T) * hedgingerror ./ V0;

% monthly hedge
% Initial portfolio
S_MC = S0 * ones(M,1); 
B = (V0-delta_bin(1,1)*S0)*ones(M,1);
delta_hedge_old = delta_bin(1,1)*ones(M,1);
delta_hedge_new = delta_hedge_old;
dt_monthly = 20*dt;
for n = 20:20:N
    % at time t_n
    % obtain phi_n for each path
    sample = randn(M,1);
    % calculate S(t_{n+1}), see (8) in assignment 2
    S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt_monthly + sigma.* sample .* sqrt(dt_monthly));
    delta_hedge_new = interpDelta(delta_bin(1:n+1, n+1), S(1:n+1, n+1), S_MC);
    B = B.*exp(r * dt_monthly) + (delta_hedge_old - delta_hedge_new).* S_MC;
    delta_hedge_old = delta_hedge_new;
end
% note that last entry of 20:20:N is 240
% to liquidate the porfolio, we need to use dt_10 = 10*dt

% at t_N, N = 250
dt_10 = 10*dt;
sample = randn(M,1);
% calculate S(t_{n+1}), see (8) in assignment 2
S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt_10 + sigma.* sample .* sqrt(dt_10));
W = max(S_MC -K, 0) + max(K- S_MC, 0);
hedgingerror = -W+ delta_hedge_new .* S_MC + B.*exp(r*dt_10);

PL_monthly = exp(-r*T) * hedgingerror ./ V0;

%% Data
PL = [PL_nohedging, PL_daily, PL_weekly, PL_monthly];

%% compute and report the performance measures of different rebalancing times
beta = 0.95;

PLperformance_table = zeros(4,4);

for i = 1:4
    PLperformance_table(1,i) = mean(PL(:,i));
    PLperformance_table(2,i) = std(PL(:,i));
    [var,cvar] = dVaRCVaR(PL(:,i),beta);
    PLperformance_table(3,i) = var;
    PLperformance_table(4,i) = cvar;
end

PLperformance_table = array2table(PLperformance_table', 'VariableNames',...
    {'Mean','Standard deviation','VaR(95%)','CVaR(95%)'}, 'RowNames',...
    {'no hedging P&L', 'daily hedging P&L','weekly hedging P&L', 'monthly hedging P&L'})

%% Discussion
% Note that VaR describes the predicted minimum profit (maximum loss) with
% a specified probability conidence level over a certain period of time.
% That means that we want to have less VaR (in absolute value since our VaR
% is applied to P&L, not on loss). Similarly, CVaR is the is the average
% P&L, given P&L is less than VaR. In other words, CVaR si measuring the
% average amount of money lost given that we are in the worst (1-beta) = 5%
% scenarios. So we also want CVaR be small in absolute value.
%
% For mean, we want the mean be close to 0 since P&L being 0 means a
% perfect hedge. If the mean is close to 0, we want the standard deviation
% be small so that P&L along each path is close to 0.
%
% Rebalancing frequency: daily > weekly > monthly > no
% 
% Based on the table above, we see that the means of no hedging, daily,
% weekly and monthly hedging are all close to 0, with the means of daily
% and weekly hedging particularly close to 0, being one order better than
% no hedging and monthly hedging. In fact the means (in absolute value) is
% a decreasing function of rebalancing frequency. We also see that the
% standard deviation is a decreasing function of rebalancing frequency.
% Similarly, VaR and CVaR (in absolute value) are decreasing functions of
% rebalancing frequency. Thus, the more frequent you rebalance, the better
% hedging preformance you get. In particular, among the rebalancing times
% we consider, daily hedging has the best performance.