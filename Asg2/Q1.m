%% paremeters initialization 

sigma = 0.2;
r = 0.03;
T = 1;
K = 10;
S0 = 10;
D0 = 0.5;
t_d = T/3;

dt = 0.05; 

%% Q1a 

testlength = 9;

dtlist = (0.05 ./ 2.^(0:testlength))';

% convergence test 

% no dividends, call option
rho = 0;
D0 = 0;
opttype = 0;

convertestcall = convergencetest(sigma, r, T, K, S0, opttype, rho, D0, t_d, dt, testlength);
convertestcall = array2table(convertestcall, ...
    'VariableNames',{'Delta t','Value','Change', 'Ratio'})

% no dividends, put option 
rho = 0;
D0 = 0;
opttype = 1;

convertestput = convergencetest(sigma, r, T, K, S0, opttype, rho, D0, t_d, dt, testlength);
convertestput = array2table(convertestput, ...
    'VariableNames',{'Delta t','Value','Change', 'Ratio'})

[exactCall,exactPut] = blsprice(S0, K, r, T, sigma)

% We see that our binomial pricing value converges to the exact solution.
% Moreover, the ratio is roughly 2, indicating a linear convergence rate.

%% Q1b
sigma = 0.2;
r = 0.03;
T = 1;
K = 10;
S0 = 10;
D0 = 0.5;
t_d = T/3;

dt = 0.005;

% we add rho = 10% to show the trend of call and put values with respect to
% rho more clearly.
rholist = [0, 0.02, 0.04, 0.08,0.1];

dividendresult = zeros(length(rholist),2);

for i = 1:length(rholist)
    % call
    dividendresult(i,1) = mybin_div(sigma, r, T, K, S0,dt, 0, rholist(i), D0, t_d);
    % put
    dividendresult(i,2) = mybin_div(sigma, r, T, K, S0,dt, 1, rholist(i), D0, t_d);
end

dividendresult = array2table(dividendresult', 'VariableNames',...
    {'rho=0%', 'rho=2%','rho=4%','rho=8%','rho=10%'}, 'RowNames',{'Call','Put'})

% We see that call values decrease as dividend yield rho increases. Also,
% we see that put values increase as dividend yield rho decreases. It can
% intutively explained. A higher dividend yield rho means that the
% ex-dividend stock value will be smaller and hence the possible payoff
% (thinking in terms of expectation) for call option will decrease while
% the possible payoff for put option (thinking in terms of expectation)
% will increase.


function value = mybin_div(sigma, r, T, K, S0,dt, opttype, rho, D0, t_d)
% mybin_div return put/call options value with a discrete dividend
%   More explanation on the formula for dividend can be found in the
%   assignment.
% 
% Input: 
% sigma: volatility of the underlying 
% r: interest rate
% T: time of expiry
% K: strike price
% dt: size of timestep 
% opttype: option type, 0 for call, otherwise for put 
% rho: constant dividend rate
% D0: constant dividend floor 
% t_d: dividend payment time

% up and down ratio in bin. model
u = exp(sigma * sqrt(dt));
d = exp(-sigma * sqrt(dt));
% risk neutral probability of having an up
q = (exp(r * dt) - d) / (u-d);

% find the index of the closest timestep that is larger than t_d
err = (0:dt:T) - t_d ;
idx = find(err>= 0, 1, 'first'); 
% so (idx - 1)th timestep (i.e. (idx-1)*dt) approximates t^+, note that we
% count the timestep from 0 to N

N = T / dt; % assume it is an integer 

% vectorized approach, find payoff at final time T, denoted by W 

% first: the stock values at final time T
%           values are arranged in desceding order 
S = S0*d.^([N:-1:0]') .* u.^([0:N]');

% second: distinguish the case between call and put
if opttype == 0
    W = max(S - K,0);
else 
    W = max(K - S,0);
end

% Backward iteration
for i = N:-1:1
    S = S(2:i+1,1) / u; % asset prices at (i-1)th timestep (i.e. (i-1)*dt)
                        % we count timestep from 0 to N
    if i == idx % when we are at the timestep approximation of t^+
        
        % note that the values of S right now actually refers to S(t^-)
        % since no deduction on S has been made yet and we know 
        % S(t^+) = S(t^-) - D. Thus, based on dividend formula:
        div_value = max(rho * S, D0);
        
        % obtain option values at (i-1)th timestep, by risk neutral valuation
        W = exp(-r *dt) *(q* W(2:i+1) + (1-q)* W(1:i));
        
        % compute option values at t^- by interpolation 
        W = dividend(W, S, div_value);
    else
        % obtain option values at (i-1)th timestep, by risk neutral valuation
        W = exp(-r *dt) *(q* W(2:i+1) + (1-q)* W(1:i));
    end
end

value = W(1);
end

function W_out = dividend( W_in, S, div_value)
% W_in: value of option at t^+
% S : asset prices (at t^-)
% div_value: discrete dollar dividend
%
% W_out: option value at t^-
%
% assume American constraint applied in caller

S_min = min(S);
S_ex = S - div_value; % ex dividend stock value
S_ex = max( S_ex, S_min); % make sure that
                            % dividend payment does
                            % not cause S < S_min
W_out = interp1( S, W_in, S_ex);
end

function testtable = convergencetest(sigma, r, T, K, S0, opttype, rho, D0, t_d, dt, testlength)
% convergencetest return a table that contains results of convergence test
% as described by Table 2 in assignment 2 
% 
% Similar input parameters as mybin_div, with the additional parameters dt
% and testlength.
% dt: initial size of testing timestep 
% testlength: (the number of size of timesteps) - 1, note that each subsequent
% size of timestep is obtained by size of previous timestep divided by 2

dtlist = (dt ./ 2.^(0:testlength))';

testtable = zeros(testlength + 1, 4);

testtable(:,1) = dtlist;
for i = 1:testlength + 1
    testtable(i,2) = mybin_div(sigma, r, T, K, S0,dtlist(i), opttype, rho, D0, t_d);
end

testtable(2:end,3) = testtable(2:end,2) - testtable(1:end-1,2);
testtable(3:end,4) = testtable(2:end-1,3)./ testtable(3:end,3);
end

