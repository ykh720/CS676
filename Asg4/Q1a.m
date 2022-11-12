%% Q1a

clearvars
close all
% parameters initialization
S0 = 76;
T = 1;
Deltatau = T/25;
alpha = 15;
r = 0.03;
K = S0;

% check S 
S = [ 0:0.1*S0:0.4*S0, ...
    0.45*S0:0.05*S0:0.8*S0, 0.82*S0 :0.02*S0 :0.9 *S0,...
    0.91*S0 :0.01*S0 :1.1*S0, 1.12*S0 :0.02*S0 :1.2*S0, ...
    1.25*S0 :0.05*S0 :1.6*S0, 1.7*S0 :0.1*S0 :2*S0,...
    2.2*S0, 2.4*S0, 2.8*S0, 3.6*S0, 5*S0, 7.5*S0, 10*S0];
S = S';
Scell = {};
Scell{1} = S;

Vold = max(K-S,0).^3;
figure(1)
plot(S, Vold)
title('payoff graph, power put')
xlabel('asset value S')
ylabel('payoff')

Deltataulist = zeros(5,1);
Deltataulist(1) = Deltatau;
% 4 refinements
for i = 1:4
    Sold = Scell{i};
    Snew = [Sold;(Sold(1:end-1) + Sold(2:end)) / 2 ];
    Snew = sort(Snew);
    Scell{i+1} = Snew;
    Deltataulist(i+1) = Deltatau/(2^i);
end

% make convergence table
Fimplicittable = zeros(length(Deltataulist),5);
CNtable = zeros(length(Deltataulist),5);
CNRantable = zeros(length(Deltataulist),5);

for i = 1:5
    % obtain option value at t=0, S = S0 using the corresponding grid and
    % time step
    [V0FL, V0CN,V0CNRan] = powerputpricer(Scell{i}, Deltataulist(i), S0, T, K, r);
    
    % nodes number 
    Fimplicittable(i,1) = size(Scell{i},1);
    CNtable(i,1) = size(Scell{i},1);
    CNRantable(i,1) = size(Scell{i},1);
    
    % timesteps number
    Fimplicittable(i,2) = 25 * 2^(i-1);
    CNtable(i,2) = 25 * 2^(i-1);
    CNRantable(i,2) = 25 * 2^(i-1);
    
    % option Value
    Fimplicittable(i,3) = V0FL;
    CNtable(i,3) = V0CN;
    CNRantable(i,3) = V0CNRan;
    
    % change in option value 
    if i > 1
        Fimplicittable(i,4) = Fimplicittable(i,3) - Fimplicittable(i-1,3);
        CNtable(i,4) = CNtable(i,3) - CNtable(i-1,3);
        CNRantable(i,4) = CNRantable(i,3) - CNRantable(i-1,3);
    end
    
    % Ratio
    if i > 2 
        Fimplicittable(i,5) = Fimplicittable(i-1,4)/Fimplicittable(i,4);
        CNtable(i,5) = CNtable(i-1,4)/CNtable(i,4);
        CNRantable(i,5) = CNRantable(i-1,4)/CNRantable(i,4);
    end 
end

Fimplicittable = array2table(Fimplicittable, 'VariableNames', ...
    {'Nodes','timesteps','Value','Change','Ratio'})
CNtable = array2table(CNtable, 'VariableNames', ...
    {'Nodes','timesteps','Value','Change','Ratio'})
CNRantable = array2table(CNRantable, 'VariableNames', ...
    {'Nodes','timesteps','Value','Change','Ratio'})

%% Q1a Discussion
% Note that fully implicit method has convergence rate O(\Delta tau, (\Delta S)^2)
% and hence it does not satisfy (4) in our question. Thus we don't expect
% the ratio for the convergence table of fully implicit method to be 4. In
% fact, we observe that such ratio is around 3, reflecting that the
% convergence rate is slightly better than first order (since we have
% O((\Delta S)^2)) but not as good as second order (since we have O(\Delta
% tau)). 
% 
% Note that Crank-Nicolson method has convergence rate  O((\Delta tau)^2,
% (\Delta S)^2). Hence, we expect the ratio for the convergence table of CN
% method to be 4. We observe that such ratio is indeed roughly 4,
% consistent with the theory regarding the rate of convergence. 
% Note that generally non-smooth payoff will cause oscillations in
% solution, resulting in slow convergence of CN method. However, in our
% case, the payoff function of the power put is quite smooth since we are
% taking cubic power. The smoothness of the payoff can also be seen in
% figure 1. At least, the payoff is continuous and the fact that K = S0 and
% we have a node at S0 = K means that no smoothing is required and we
% expect the usual second order convergence.
% 
% For CN-Rannacher, the purpose of it is to do some smoothing so that
% second order convergence is attained even for nonsmooth payoff. But in
% our case, the payoff is quite smooth. But it does not prohibit
% CN-Rannacher to achieve second order convergence, as demostrated in the
% Ratio column of its convergence table being around 4.

function [V0fullyimplicit, V0CN,V0CNRan] = powerputpricer(S, Deltatau, S0, T, K, r)
% output option value at t =0, S =S0, using fully implicit, Crank Nicolson,
% and CN-Rannacher

%%% Fully implicit
%%% only compute initial option value
%%% can easily modified to obtain options values at each grid points
m = length(S);
Vold = zeros(m,1); 
Vnew = Vold;
N = T / Deltatau;
% initialize payoff 
Vold = max(K-S,0).^3;

% preprocessing step for simple BS equation
% so that we have positive coefficient discretization
% these will also be used for CN and CN-Rannacher
Smid = S(2:end-1);
Sleft = S(1:end-2);
Sright = S(3:end);
sigma = sigmalvf(Smid);
alphacentral = sigma.^2 .* Smid .^2 ./ ((Smid - Sleft).* (Sright - Sleft)) - ...
    r* Smid ./(Sright - Sleft);
betacentral = sigma.^2 .* Smid .^2 ./ ((Sright - Smid).* (Sright - Sleft)) + ...
    r* Smid ./(Sright - Sleft);

alphaups = sigma.^2 .* Smid .^2 ./ ((Smid - Sleft).* (Sright - Sleft));
betaups = sigma.^2 .* Smid .^2 ./ ((Sright - Smid).* (Sright - Sleft)) + ...
    r* Smid ./(Sright - Smid);

idx1 = alphacentral <0;
idx2 = betacentral <0;
idx = idx1 | idx2;

alphavec = alphacentral;
betavec = betacentral;
alphavec(idx) = alphaups(idx);
betavec(idx) = betaups(idx);

% build the matrix M by spdiags
Bin = Deltatau* [[[-alphavec;0];0], [[r; (alphavec + betavec +r)];0],[[0;0];-betavec]];
M = spdiags(Bin, -1:1, m, m);
% spy(M) % see the nonzero entries of the matrix M
% Note that our M is time independent
I = speye(m);
[L,U,P] = lu(M+I); 
% increase performance since the LU factorization is only computed once and
% being reused in each linear system solving
% to solve (M+I)x = b, y =L\ (P*b), x = U\y

for n = 1 : N  
    y = L\(P*Vold);
    Vnew = U\y;
    Vold = Vnew;
end
S0id = S0 == S;
V0fullyimplicit = Vnew(S0id);


%%% Crank-Nicolson
%%% only compute initial option value
%%% can easily modified to obtain options values at each grid points
Vold = zeros(m,1); 
Vnew = Vold;
% initialize payoff 
Vold = max(K-S,0).^3;

Bin = Deltatau/2* [[[-alphavec;0];0], [[r; (alphavec + betavec +r)];0],[[0;0];-betavec]];
Mhat = spdiags(Bin, -1:1, m, m);
% Note that our Mhat is time independent
I = speye(m);
[L1,U1,P1] = lu(Mhat+I); 


for n = 1 : N  
    y = L1\(P1*(I-Mhat)*Vold);
    Vnew = U1\y;
    Vold = Vnew;
end
V0CN = Vnew(S0id);

%%% CN-Rannacher
Vold = zeros(m,1); 
Vnew = Vold;
% initialize payoff 
Vold = max(K-S,0).^3;

% two steps of fully implicit 
for n = 1 : 2  
    y = L\(P*Vold);
    Vnew = U\y;
    Vold = Vnew;
end

% use CN after that
for n = 1 : N -2
    y = L1\(P1*(I-Mhat)*Vold);
    Vnew = U1\y;
    Vold = Vnew;
end
V0CNRan = Vnew(S0id);

end

function out = sigmalvf(S)
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
    alpha = 15;
    S0 = 76;
    out = alpha./ (max(S0 + (S- S0)./ (2*S0), 10^(-8)));
end