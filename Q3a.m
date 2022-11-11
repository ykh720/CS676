%% Parameters initialization 

K = 102;
B = 100;
T = 0.5;
x = 15;
sigma = 0.2;
r = 0.03;
S0 = 100:2:120;

% fair value of down-and-out cash-or-nothing put option 
V = do_cashornth_put(sigma, r, T, K, S0, x, B)

% why the first entry is actually negative ??
% is it because of floating point error?, maybe it should be small very
% small positive value, but after rounding it becomes negative

% plot the down-and-out cash-or-nothing put option for S=100:2:120
h= figure;
plot(S0,V)
xlabel('initial stock value: S0')
ylabel('option value: V')
saveas(h, 'Q3a','epsc') % stored as EPS, there is no loss in quality 


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
