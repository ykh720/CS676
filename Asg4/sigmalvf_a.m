function sigma_a = sigmalvf_a(S)
%%  local volatility by \alpha ./ sqrt(S)
% 
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
%
% inputs: S, stock A price
alpha = 0.85;
sigma_a = alpha./ sqrt(S);
end