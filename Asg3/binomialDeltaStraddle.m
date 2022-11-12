% Q3a
function [V0,S, delta] = binomialDeltaStraddle(S0,r,sigma, T,N,K)
% V0 option price at time 0 
% S, matrix of underlying prices, the i-th column represent t_{i+1}, only
% store until t_{N-1}
% delta, matrix of delta, same storage as S 

dt = T/N;
% up and down ratio in bin. model
u = exp(sigma * sqrt(dt));
d = exp(-sigma * sqrt(dt));
% risk neutral probability of having an up
q = (exp(r * dt) - d) / (u-d);

S = zeros(N,N);
delta = zeros(N,N);

% the stock values at final time T
%           values are arranged in ascending (from top to bottom) order 
Svec = S0*d.^([N:-1:0]') .* u.^([0:N]');

% final payoff 
W = max(Svec - K,0) + max(K - Svec,0);

% fill in S and delta, column by column 
for i = N:-1:1
    Svec = Svec(2:i+1) ./ u;
    S(1:i,i) = Svec;
    delta(1:i,i) = (W(2:end) - W(1:end-1)) ./ ((u-d)*Svec);
    % obtain option values at (i-1)th timestep, by risk neutral valuation
    W = exp(-r *dt) *(q* W(2:i+1) + (1-q)* W(1:i));
end

V0 =W;
end

