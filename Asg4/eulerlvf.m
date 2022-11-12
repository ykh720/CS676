% Q6a 
function Snew = eulerlvf(Sold, T, N, opt, r)
% Sold is the stock price at time 0
% Snew will be the stock price at time T, generated by Euler method
%
% opt = 1 for stock A, else for stock B
% N: number of timesteps taken, integer, so dt = T / N
%
% assume column vector as input
m= length(Sold);
Snew = zeros(m,1);
dt = T/N;
sqrt_dt = sqrt(dt);

if opt == 1 
    for i = 1: N
        phi = randn(m,1);
        Snew = Sold + Sold .*(r * dt + sigmalvf_a(Sold)* sqrt_dt .*phi);
        Sold = Snew;
    end
else 
    for i = 1: N
        phi = randn(m,1);
        Snew = Sold + Sold .*(r * dt + sigmalvf_b(Sold)* sqrt_dt .*phi);
        Sold = Snew;
    end 
end

end
