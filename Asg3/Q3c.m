%Q3c

sigma = 0.2;
r = 0.03;
mu = 0.15;
T = 1;
S0 = 10;
K = 0.95 * S0;
N = 250;

% query points
S_query = linspace(0.8*S0, 1.4*S0, 100);

[V0,S, delta] = binomialDeltaStraddle(S0,r,sigma, T,N,K);

n = 0.8*N;
delta_n = delta(1:n+1,n+1);
S_n = S(1:n+1, n+1);

delta_query = interpDelta(delta_n,S_n, S_query);

h = figure;
plot(S_query, delta_query)
title('hedge position delta against S')
xlabel('S')
ylabel('delta')
saveas(h, 'Q3c','epsc')