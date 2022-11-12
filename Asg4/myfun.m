%% Q3c

function [F,J] = myfun(x)
% marketvalue: V_0^mkt, a column vector 

% parameters and data preparation
S0 =100;
r = 0.03;
Klist = [0.9, 0.95, 1, 1.05, 1.1]* S0;
Tlist = [0.425, 0.695, 1];
N = 100; % what should be the correct number of steps?
% can change to N = 100

S = [ 0:0.1*S0:0.4*S0, ...
    0.45*S0:0.05*S0:0.8*S0, 0.82*S0 :0.02*S0 :0.9 *S0,...
    0.91*S0 :0.01*S0 :1.1*S0, 1.12*S0 :0.02*S0 :1.2*S0, ...
    1.25*S0 :0.05*S0 :1.6*S0, 1.7*S0 :0.1*S0 :2*S0,...
    2.2*S0, 2.4*S0, 2.8*S0, 3.6*S0, 5*S0, 7.5*S0, 10*S0];
S = S';

% can do the refinements following 
Sold = S;

for i = 1:2
    Snew = [Sold;(Sold(1:end-1) + Sold(2:end)) / 2 ];
    Snew = sort(Snew);
    Sold = Snew;
end

S = Snew;


% IV table
table2 = zeros(3,5);
table2(1,:) = [0.155, 0.138, 0.125, 0.109, 0.103];
table2(2,:) = [0.157, 0.144, 0.133, 0.118, 0.104];
table2(3,:) = [0.159, 0.149, 0.137, 0.127, 0.113];

market = zeros(3,5);

for j = 1:5
    for i = 1:3
        [Call , ~] = blsprice(S0, Klist(j), r, Tlist(i), table2(i,j));
        market(i,j) = Call;
    end
end

market = reshape(market, numel(market), []);

F = zeros(length(Tlist),length(Klist));  
% Objective function values at x, stored as in table 2

S0id = S == S0;

for i = 1:length(Tlist)
    T = Tlist(i);
    Deltatau = T/N;
    for j = 1 : length(Klist)
        K = Klist(j);
        V0CNRan = callpricer_CNR(S, Deltatau, x, T, K, r);
        % only consider CNR, note that other methods can be obtained in
        % callpricer.m
        value = V0CNRan(S0id);
        F(i,j) = value;
    end
end

% Note that F now only contains value of options depending on x, not
% subtracted by the the market values yet

% reshape F, column by column
F = reshape(F, numel(F),[]);


% Note that we don't need to consider the constant market values when doing
% Jacobian. Hence we will subtract F by marketvalue later.
if nargout > 1 % Two output arguments
    change = sqrt(eps);
    Fnew = zeros(length(Tlist),length(Klist));
    J = zeros(numel(Fnew), length(x)); % Jacobian of the function evaluated at x
    
    for k = 1 : length(x)
        Fnew = zeros(length(Tlist),length(Klist));
        % calculate the k-th column of J
        dir = zeros(length(x), 1); % direction of derivative
        dir(k) = 1;
        for i = 1:length(Tlist)
            T = Tlist(i);
            Deltatau = T/N;
            for j = 1 : length(Klist)
                K = Klist(j);
                V0CNRan = callpricer_CNR(S, Deltatau, x+ change*dir, T, K, r);
                % only consider CNR, note that other methods can be obtained in
                % callpricer.m
                Fnew(i,j) = V0CNRan(S0id);
            end
        end
        Fnew = reshape(Fnew, numel(Fnew), []);
        J(:,k) = (Fnew - F) / change;
    end
        
end

F = F - market;

end
