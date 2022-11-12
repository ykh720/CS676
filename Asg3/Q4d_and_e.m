%% Q4d_and_e

close all 
clearvars

rng('default')

% parameter initialization

M = 80000; 
sigma = 0.2;
r = 0;
mu = 0.15;
T = 1;
S0 = 100;
K = 105;
Nlist = [100,200,400,800];

% [V0,S, delta_bin] = binomialDeltaStraddle(S0,r,sigma, T,Nlist(end),K);

[exactcall, exactput] = blsprice(S0, K, r, T, sigma)
PL = zeros(M,4);
for i = 1:length(Nlist)
    N = Nlist(i);
    % disp(N)
    % initial positions
    dt = T/N;
    % disp(dt)
    S_MC = S0 * ones(M,1); 
    B = (exactput + S0)*ones(M,1);
    delta_hedge_old = -ones(M,1);
    delta_hedge_new = delta_hedge_old;
    
    % it seems like we also rebalance at the end t_N so we use 1:N instead of
    % 1:N-1
    for n = 1:N
        % logical values from t_{n-1}
        idx1 = S_MC < K; % delta_t = 1 if S_t >= K, although equality should not happen 
        idx2 = ~idx1; % S_MC >= K

        % at time t_n
        % obtain phi_n for each path
        sample = randn(M,1);

        % calculate S(t_{n+1}), see (8) in assignment 2, using mu
        S_MC = S_MC.*exp((mu - 1/2 * sigma^2) * dt + sigma.* sample .* sqrt(dt));

        % update delta only when (S_n < K and S_{n+1} >= K) or (S_n >= K and
        % S_{n+1} < K), the former case update delta to 0, while the latter
        % update it to -1

        idx3 = idx1 & (S_MC >= K);
        idx4 = idx2 & (S_MC < K);

        delta_hedge_new(idx3) = 0; % 0 or 1???
        delta_hedge_new(idx4) = -1;

        B = B.*exp(r * dt) + (delta_hedge_old - delta_hedge_new).* S_MC;
        delta_hedge_old = delta_hedge_new;
    end
    W = max(K- S_MC,0);

    PL(:,i) = (-W + delta_hedge_new .* S_MC + B)./ (S0 + exactput);
end

beta = 0.95;

PLperformance_table = zeros(4,4);

for i = 1:4
    PLperformance_table(1,i) = mean(PL(:,i));
    PLperformance_table(2,i) = std(PL(:,i));
    [var,cvar] = dVaRCVaR(PL(:,i),beta);
    PLperformance_table(3,i) = var;
    PLperformance_table(4,i) = cvar;
end

PLperformance_table = array2table(PLperformance_table, 'RowNames',...
    {'Mean','Standard deviation','VaR(95%)','CVaR(95%)'}, 'VariableNames',...
    {'100', '200','400', '800'})

% probability density plot
h = figure(1);
[f,xi] = ksdensity(PL(:,4));
plot(xi,f)
title('probability density plot for 800 rebalancing times')
xlabel('relative P&L')
ylabel('pdf')
saveas(h,'Q4fig1','epsc')

% histogram with normal fit
g = figure(2);
histfit(PL(:,4), 50)
title('histogram for 800 rebalancing times with normal fit')
saveas(g,'Q4fig2','epsc')

%% Discussion
% What do you observe about the mean and variance of the hedging error?
% I observe that the mean is actually less than 0, although being quite
% close to 0 (around -0.0005). The standard deviation is around 0.053,
% which is quite low. We see that for there are no substantial difference
% bewteen different rebalancing times. They all have similar (at least same
% order of magnitude) mean, standard deviation, VaR and CVaR. 

%% Discussion, part e
% I have also tried to use some very high rebalancing times like 30000. But
% the mean of P&L is still around -0.0052306, which is roughly the same for
% rebalancing times 100,200,400,800. It seems to me that for each
% transcation, we incur some small losses since we are selling the stock
% with price < K and buying the stock with price > K, owing to the time
% discretization error so that we cannot exactly buy and sell the stock at
% price K. If we refine our timestepping, i.e. taking larger N, then the
% loss for each transaction is smaller but we also have to do more
% transactions. So that the total loss stays around the same. That might
% possibly explains why even we use a much higher rebalancing time, there
% is still no substantial improvment. This discussion also addresses part
% e.