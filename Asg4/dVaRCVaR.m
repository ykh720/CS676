%% ask, the input is profit or loss? 
% after reading lecture carefully, i think it is profit

function [var,cvar] = dVaRCVaR(PL, beta)
% PL: discrete P&L distribution with M independent samples
% assume PL is a column vector
% assume PL is profit, i.e. positive value for profit
% beta could be 95% 

M = length(PL);
PL_sorted = sort(PL); % sorts in ascending order.

idx = floor((1-beta)*M);
var = PL_sorted(idx);

cvar = mean(PL_sorted(1:idx));

end

