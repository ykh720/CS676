%% Q2b
% use CNRall, a version of CNR in Q2a that outputs the whole options value
% associated to the grid at t = 0

clearvars
close all
%%  parameters initialization
S0 = 76;
T = 1;
Deltatau = T/25; % initial time step
dnorm = 0.1;
alpha = 15;
r = 0.03;
K = S0;

S = [ 0:0.1*S0:0.4*S0, ...
    0.45*S0:0.05*S0:0.8*S0, 0.82*S0 :0.02*S0 :0.9 *S0,...
    0.91*S0 :0.01*S0 :1.1*S0, 1.12*S0 :0.02*S0 :1.2*S0, ...
    1.25*S0 :0.05*S0 :1.6*S0, 1.7*S0 :0.1*S0 :2*S0,...
    2.2*S0, 2.4*S0, 2.8*S0, 3.6*S0, 5*S0, 7.5*S0, 10*S0];
S = S';

%% Now do refinements
Scell = {};
Scell{1} = S;
% list of initial timestep
Deltataulist = zeros(5,1);
Deltataulist(1) = Deltatau;
dnormlist = zeros(5,1);
dnormlist(1) = dnorm;
% 4 refinements
for i = 1:4
    Sold = Scell{i};
    Snew = [Sold;(Sold(1:end-1) + Sold(2:end)) / 2 ];
    Snew = sort(Snew);
    Scell{i+1} = Snew;
    Deltataulist(i+1) = Deltatau/(4^i);
    dnormlist(i+1) = dnorm/(2^i);
end
% use finest grid 
S = Scell{end};
Deltatau = Deltataulist(end);
dnorm = dnormlist(end);

%% load data for european, obtain data for American, plot graphs

filename = 'Q1b.mat';
myVars = {'Vat0','deltaat0','gammaat0'};

% load data from Q1b 
Euro = load(filename, myVars{:});

[V0, steptotal, counter] = CNRall(S, Deltatau, S0, T, K, r,dnorm,1);

Vat0 = V0;

% Central difference for delta and gamma at t = 0

idx = (S <= 1.5*S0) & (S >= 0.5 * S0);
idx = find(idx);
index1 = idx(1);
index2 = idx(end);

deltaxright = S(index1 + 1: index2 + 1) - S(index1:index2);
deltaxleft = S(index1:index2) - S(index1 - 1:index2-1);

deltafright0 = Vat0(index1 + 1: index2 + 1) - Vat0(index1:index2);
deltafleft0 = Vat0(index1:index2) - Vat0(index1 - 1:index2-1);
gammaat0 = (deltafright0 ./ deltaxright - deltafleft0./deltaxleft) ...
    ./ ((deltaxright + deltaxleft)/2);
% not sure about this, I think this might improve the order of delta a bit 
deltaat0 = (Vat0(index1+1: index2+1)-Vat0(index1-1:index2-1))...
    ./(deltaxright + deltaxleft) - gammaat0 /2 .*(deltaxright-deltaxleft);

figure(1)
plot(S(idx), Euro.Vat0(idx))
hold on 
plot(S(idx),Vat0(idx))
title('Option value')
xlabel('asset value S')
legend('Euro','American')
hold off

figure(2)
plot(S(idx), Euro.deltaat0)
hold on 
plot(S(idx), deltaat0)
title('Delta')
xlabel('asset value S')
legend('Euro','American')
hold off

figure(3)
plot(S(idx), Euro.gammaat0)
hold on 
plot(S(idx), gammaat0)
title('Gamma')
xlabel('asset value S')
legend('Euro','American')
hold off

figure(4) 
plot(S(idx), (Vat0(idx)- Euro.Vat0(idx)) )
title('American - European')
xlabel('asset value S')

figure(5)
plot(S(idx), (Vat0(idx)- Euro.Vat0(idx))./Euro.Vat0(idx) )
title('American - European, relative change')
xlabel('asset value S')

%% Discussion
% Generally, we see that the European options value is smaller than or
% equal to  that of American options. This is because American options can
% exercise early and hence its value must be >= European options. However,
% such difference is very small in our case. This is because r = 0.03 is
% quite small and the interval [0.5 *S0, 1.5*S0] that we consider is not
% very at the money for the power put options, so that the early exercise
% value is not very large. But note that the as S decreases (and so payoff
% increases), such early exercise value is larger and a larger difference
% is thus shown on the left of the option value plot.
% 
% Note that if we set r = 0.5, then a much bigger difference will be seen.
% If we set r = 0, then no difference will be seen.
