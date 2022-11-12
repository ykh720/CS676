function sigma_b = sigmalvf_b(S)
%%  local volatility by \alpha ./ sqrt(S)
% 
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
%
% inputs: S, stock B price
alpha = 0.9;
sigma_b = alpha./ sqrt(S);
end