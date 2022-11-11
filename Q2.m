% This is the latest script for Q2, more precise, less comprehensive.

%% paremeters initialization 

sigma = 0.2;
r = 0.03;
T = 1;
K = 10;
S0 = 10;
dt = 0.05; 

%% Q2a 

testlength = 9;

% Exact solution 
[exactCall,exactPut] = blsprice(S0, K, r, T, sigma)

% put option
opttype = 1;

convertestput = convergencetest_drift(sigma, r, T, K, S0, opttype, dt, testlength);
convertestput = array2table(convertestput, ...
    'VariableNames',{'Delta t','Value','Change', 'Ratio'})

% We see that the put values converge to the exact solutions as Delta t
% goes to 0.

% We observe that the ratio in the above table are not constant and
% flutuate a lot. This means that the convergenece rate is neither linear
% nor quadartic. By p.11 in Lec 5, this may indicate that the strike price
% is not at a binomial mode.

%% Q2b smoothed payoff
% call option 
opttype = 0;
smoothtestcall = convergencetest_drift_smoothed(sigma, r, T, K, S0, opttype, dt, testlength);
smoothtestcall = array2table(smoothtestcall, ...
    'VariableNames',{'Delta t','Value','Change', 'Ratio'})

% put option 
opttype = 1;
smoothtestput = convergencetest_drift_smoothed(sigma, r, T, K, S0, opttype, dt, testlength);
smoothtestput = array2table(smoothtestput, ...
    'VariableNames',{'Delta t','Value','Change', 'Ratio'})

function value = driftlat(sigma, r, T, K, S0,dt, opttype)
% driftlat return put/call options value using drifting lattice
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


term = (r- sigma^2/2 ) * dt ;
% up and down ratio in drifting lattice
u = exp(sigma * sqrt(dt) + term);
d = exp(-sigma * sqrt(dt) + term);
% probability of going up
q = 1/2;

N = T / dt; % assume it is an integer 

% vectorized approach, find payoff at final time T 
% first: the stock values at final time T
%           values are arranged in desceding order 
W = S0*d.^([N:-1:0]') .* u.^([0:N]');

% second: distinguish the case between call and put
if opttype == 0
    W = max(W - K,0);
else 
    W = max(K - W,0);
end

% Backward iteration
for i = N:-1:1
    W = exp(-r *dt) *(q* W(2:i+1) + (1-q)* W(1:i));
end

value = W(1);
end

function value = driftlat_smoothed(sigma, r, T, K, S0,dt, opttype)
% similar to driftlat, except that we implement smoothing payoff

term = (r- sigma^2/2 ) * dt ;
u = exp(sigma * sqrt(dt) + term);
d = exp(-sigma * sqrt(dt) + term);
q = 1/2;

N = T / dt; % assume it is an integer 

% vectorized approach, find payoff at final time T 
% first: the stock values at final time T
%           values are arranged in desceding order
W = S0*d.^([N:-1:0]') .* u.^([0:N]');

% some useful constants to simplify the expression
up = exp(sigma * sqrt(dt));
down = exp(-sigma * sqrt(dt));

% second: distinguish the case between call and put
% smoothed payoff, cf (5.49), (5.50) in the pdf course notes
if opttype == 0
    % call case:
    idx = (W *down <= K) & (K <= W*up);
    idx1 = W* down > K;
    idx2 = W *up < K;
    W(idx2) = 0;
    W(idx1) = W(idx1).*(up -down)./ (2 * sigma * sqrt(dt)) - K;
    W(idx) = (1 / (2* sigma * sqrt(dt))) .*(W(idx) .*(up - K./W(idx))- ...
        K *(sigma * sqrt(dt) - log(K./ W(idx))));
else 
    % put case:
    idx = (W *down <= K) & (K <= W*up);
    idx1 = W* down > K;
    idx2 = W *up < K;
    W(idx1) = 0;
    W(idx2) = K - W(idx2).*(up -down)./ (2 * sigma * sqrt(dt));
    W(idx) = (1 / (2* sigma * sqrt(dt))) .*(K *(log(K./ W(idx)) + sigma * sqrt(dt)) ...
        - W(idx) .*(K./W(idx) - down));
end

% Backward iteration
for i = N:-1:1
    W = exp(-r *dt) *(q* W(2:i+1) + (1-q)* W(1:i));
end

value = W(1);
end

function testtable = convergencetest_drift(sigma, r, T, K, S0, opttype, dt, testlength)
% convergencetest_drift return a table that contains results of convergence test
% for drifting lattice (without smoothing payoff)
% 
% Similar input parameters as driftlat, with the additional parameters dt
% and testlength.
% dt: initial size of testing timestep 
% testlength: (the number of size of timesteps) - 1, note that each subsequent
% size of timestep is obtained by size of previous timestep divided by 2


dtlist = (dt ./ 2.^(0:testlength))';

testtable = zeros(testlength + 1, 4);

testtable(:,1) = dtlist;
for i = 1:testlength + 1
    testtable(i,2) = driftlat(sigma, r, T, K, S0,dtlist(i), opttype);
end

testtable(2:end,3) = testtable(2:end,2) - testtable(1:end-1,2);
testtable(3:end,4) = testtable(2:end-1,3)./ testtable(3:end,3);
end

function testtable = convergencetest_drift_smoothed(sigma, r, T, K, S0, opttype, dt, testlength)
% convergencetest_drift_smoothed return a table that contains results of convergence test
% for drifting lattice with smoothing payoff
% 
% Similar input parameters as driftlat_smoothed, with the additional parameters dt
% and testlength.
% dt: initial size of testing timestep 
% testlength: (the number of size of timesteps) - 1, note that each subsequent
% size of timestep is obtained by size of previous timestep divided by 2

dtlist = (dt ./ 2.^(0:testlength))';

testtable = zeros(testlength + 1, 4);

testtable(:,1) = dtlist;
for i = 1:testlength + 1
    testtable(i,2) = driftlat_smoothed(sigma, r, T, K, S0,dtlist(i), opttype);
end

testtable(2:end,3) = testtable(2:end,2) - testtable(1:end-1,2);
testtable(3:end,4) = testtable(2:end-1,3)./ testtable(3:end,3);
end



