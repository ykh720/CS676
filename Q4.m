close all

%% parameters initialization
% tstart = tic; & count the time spent 

dtlist = [0.02, 0.01, 0.005, 0.0025];
M = 100000; % quite large 

K = 102;
B = 100;
T = 0.5; % expiry of 6 months
x = 15;
sigma = 0.2;
r = 0.03; % use 0.03 or 0.1? (paper value: 0.1), (Table 1 in asg: 0.03)
% based on discussion on piazza, choose 0.03

% some numeric experiements show that if we change sigma to be less than r,
% we will have large error and that the modified MC actually not as good as
% usual MC. This can be seen by choosing r = 0.1, sigma = 0.1.
% I am not sure why this happens. Would it be possible that the analytic
% formula in part 3a can only be applied to the case that sigma > r? Again,
% I am not sure.

S0 = 105;

%% reproduction of figures of example 1 in the paper

rng('default') % for reproducibility 

% exact put option value
exact = do_cashornth_put(sigma, r, T, K, S0, x, B)

% V for usual MC price, V_mod for modified MC
V = zeros(length(dtlist),1);
V_mod = V;

for i = 1:length(dtlist)
    V(i) = MC_barrier(sigma, r, T, K, S0, x, B, dtlist(i), M);
end
V

for i = 1:length(dtlist)
    V_mod(i) = mod_MC_barrier(sigma, r, T, K, S0, x, B, dtlist(i), M);
end
V_mod

% tend = toc(tstart) % shows the time spent for this program

% plot the modified MC prices and exact price
h = figure(1);
plot(dtlist, V_mod','-or') 
hold on 
plot(dtlist, ones(1,length(dtlist))* exact,'-|b')
xlabel('DeltaT')
ylabel('Option price')
legend('V-Monte Carlo','Vexact')
hold off

% plot the approximation errors between the standard MC and the improved MC
g = figure(2);
plot(T./dtlist, abs(V- exact),'-+b')
hold on 
plot(T./dtlist, abs(V_mod - exact), '-or')
xlabel('Number of steps')
ylabel('Error (absolute)')
legend('Standard MC','Improve MC')
hold off

saveas(h,'Q4fig1','epsc')
saveas(g,'Q4fig2','epsc')

% more plot, also show the time discretization error
g2 = figure(3);
plot(T./dtlist, abs(V- exact),'-+b')
hold on 
plot(T./dtlist, abs(V_mod - exact), '-or')
plot(T./dtlist, sqrt(dtlist), '-')
xlabel('Number of steps')
ylabel('Error (absolute)')
legend('Standard MC','Improve MC','sqrt(dtlist)')
hold off
saveas(g2,'Q4fig3','epsc')


function V = do_cashornth_put(sigma, r, T, K, S0, x, B)
% do_cashornth_put returns the exact down-and-out cash-or-nothing put option value
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

function value = mod_MC_barrier(sigma, r, T, K, S0, x, B, dt, M)
% mod_MC_barrier returns the down-and-out cash-or-nothing put option
%   value by modified Monte Carlo, using exceedance probabiliy 
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
Sold = S0 * ones(M,1); % initial asset value, will also be used as S(t_n) 
Snew = zeros(M,1); % S(t_{n+1})

for n = 0:N-1
    % obtain phi_n for each path
    sample = randn(M,1);
    % calculate S(t_{n+1}), see (8) in assignment
    Snew = Sold.*exp((r - 1/2 * sigma^2) * dt + sigma.* sample .* sqrt(dt));
    % record the indices of path of S that is below the barrier at this timestep 
    badcheck = Snew <= B;
    idx(badcheck) = false;
    
    % exceedance probability, p_{n+1} in p.68 of the paper 
    exceedprob = exp(-2.*(B-Sold(idx)).*(B - Snew(idx))./ ...
        (sigma^2.*Sold(idx).^2 * dt));
    % unifromly distributed RV, u_n
    uniformRV = unifrnd(0,1, size(exceedprob));
    % If p_n \geq u_n, that means exceedane probability is high and we
    % regard that as S reach the barrier in the time interval (t_{n},
    % t_{n+1}). We record the indices that exceedance probability is small.
    goodcheck2 = uniformRV > exceedprob;
    % only those indices with small exceedance probability can have
    % positive payoff
    idx(idx) = goodcheck2;
    
    Sold = Snew;
end

% final check: only realization with final asset value less than K can have
% cash payout
endcheck = Snew >= K;
idx(endcheck) = false;

% compute the MC price by (9) in assignment 
V(idx) = x;
value = exp(-r*T) *mean(V);
end 

