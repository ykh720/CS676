%% Q3d

% Should compare no hedging, daily, weekly, monthly hedge effectiveness
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

%% no hedging

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

PL_nohedging_mean = mean(PL_nohedging)

%% daily hedge
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

% disp('daily hedge, mean PL')
PL_mean_daily = mean(PL_daily)

%% weekly hedge
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

PL_weekly_mean = mean(PL_weekly)

%% monthly hedge
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

PL_monthly_mean = mean(PL_monthly)

%% histogram
PL = [PL_nohedging, PL_daily, PL_weekly, PL_monthly];
xlabelslist =["no hedging P&L" "daily hedging P&L" "weekly hedging P&L" "monthly hedging P&L"];
% to have string array, need to use " instead of '

% histograms with 50 bins
h = figure(1);
for i = 1:4
    subplot(2,2,i)
    binranges= linspace(min(PL(:,i)), max(PL(:,i)), 51); % use at least 50 bins
    bincounts = histc(PL(:,i), binranges);
    bar(binranges,bincounts,'histc')
    title({'histogram for', 'relative hedging error P&L'})
    ylabel('number of occurences')
    xlabel(xlabelslist(i))
end
saveas(h, 'Q3d_1','epsc')

% histogram with normal dist fit, 50 bins

h = figure(2);
for i = 1:4
    figure(2)
    subplot(2,2,i)
    histfit(PL(:,i), 50)
    title({'histogram for relative', 'hedging error P&L, with normal dist fit'})
    ylabel('number of occurences')
    xlabel(xlabelslist(i))
end
saveas(h, 'Q3d_2','epsc')

%% Comment on your observations
% We observe that for daily, weekly and monthly hedging, the P&L
% distribution is very close to normal distrbution. This is clearly
% reflected in our plots of histogram with normal distribution fit.
%
% Moreover, compared to daily, weekly and monthly hedgings, no hedging P&L
% does not behave like a normal distribution. It is more like a skewed
% distribution, but still with mean of P&L being roughly 0. No hedging P&L
% has occurrences of extremre losses while daily, weekly and monthly
% hedging do not have. This demostrates the superiority of hedging. 
%
% Among daily, weekly and monthly hedgings, we see that daily hedging
% performs the best in the sense that it has the most occurence of nearly 0
% P&L. This can be seen by observing the number of occurences of the bin
% covering 0. Note that we have a perfect hedge when the P&L is 0. In
% particular, the number of occurences for bin covering 0 for daily hedging
% is around 7000, while that for monthly hedging is only around 5500.
%
% Note that P&L is of course not zero along each scenario path (as shown in
% histograms that we have occurence of nonzero P&L) since our rebalancing
% is not continuous so there is time discretization error. Also our
% computed delta is just an approximation of the real delta, such computed
% delta is obtained via binomial model where we rely on the option values
% as computed by the binomial model. But we know that option values from
% binomial model also have time discretization error. So this another
% source of error. 
% (This part I am unsure)
% 
