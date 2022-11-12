%% Q1b
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

% USE the finest grid for CN-Rannacher 
S = Scell{end};
Deltatau = Deltataulist(end);


% preprocessing step for simple BS equation
% so that we have positive coefficient discretization
% these will also be used for CN and CN-Rannacher
m = length(S);
N = T / Deltatau;

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

Bin = Deltatau/2* [[[-alphavec;0];0], [[r; (alphavec + betavec +r)];0],[[0;0];-betavec]];
Mhat = spdiags(Bin, -1:1, m, m);
% Note that our Mhat is time independent
I = speye(m);
[L1,U1,P1] = lu(Mhat+I); 

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
Vbeforehalf = zeros(m,1);
Vafterhalf = zeros(m,1);
Vathalf = zeros(m,1);
Vbefore0 = zeros(m,1);
for n = 1 : N -2
    y = L1\(P1*(I-Mhat)*Vold);
    Vnew = U1\y;
    Vold = Vnew;
    if n == N/2 -3   % save the data at the timestep just after t= T/2,
                     % after in terms of t
        Vbeforehalf = Vold;
    end
    if n == N/2 -1   % save the data at the timestep just before t = T/2, 
                     % before in terms of t
        Vafterhalf = Vold;
    end
    if n == N/2- 2   % save the data at the t= T/2
        Vathalf = Vold;
    end
    if n == N-2 -1   % save the data at the timestep just before t = 0
        Vbefore0 = Vold;
    end
    
end

% index for S in the range = [0.5 S0, 1.5 S0]
idx = (S <= 1.5*S0) & (S >= 0.5 * S0);
idx = find(idx);
index1 = idx(1);
index2 = idx(end);

Vat0 = Vold; % the data at t = 0

% Central difference for delta and gamma at t = T/2
deltaxright = S(index1 + 1: index2 + 1) - S(index1:index2);
deltaxleft = S(index1:index2) - S(index1 - 1:index2-1);
deltafright = Vathalf(index1 + 1: index2 + 1) - Vathalf(index1:index2);
deltafleft = Vathalf(index1:index2) - Vathalf(index1 - 1:index2-1);
gammaathalf = (deltafright ./ deltaxright - deltafleft./deltaxleft) ...
    ./ ((deltaxright + deltaxleft)/2);
% not sure about this, I think this might improve the order of delta a bit 
deltaathalf = (Vathalf(index1+1: index2+1)-Vathalf(index1-1:index2-1))...
    ./(deltaxright + deltaxleft) - gammaathalf /2 .*(deltaxright-deltaxleft);

% Central difference for delta and gamma at t = 0

deltafright0 = Vat0(index1 + 1: index2 + 1) - Vat0(index1:index2);
deltafleft0 = Vat0(index1:index2) - Vat0(index1 - 1:index2-1);
gammaat0 = (deltafright0 ./ deltaxright - deltafleft0./deltaxleft) ...
    ./ ((deltaxright + deltaxleft)/2);
% not sure about this, I think this might improve the order of delta a bit 
deltaat0 = (Vat0(index1+1: index2+1)-Vat0(index1-1:index2-1))...
    ./(deltaxright + deltaxleft) - gammaat0 /2 .*(deltaxright-deltaxleft);

figure(1)
plot(S(idx), Vathalf(idx))
hold on 
plot(S(idx),Vat0(idx))
title('Option value')
xlabel('asset value S')
legend('t= T/2','t=0')
hold off

figure(2)
plot(S(idx), deltaathalf)
hold on 
plot(S(idx), deltaat0)
title('Delta')
xlabel('asset value S')
legend('t= T/2','t=0')
hold off

figure(3)
plot(S(idx), gammaathalf)
hold on 
plot(S(idx), gammaat0)
title('Gamma')
xlabel('asset value S')
legend('t= T/2','t=0')
hold off

save Q1b.mat

%% Discussion
% For figure(1) option value plot, we observe that the option values at t =
% 0 are larger than that at t = T/2, for the same asset value S. This makes
% sense since the Greek Theta is generally negative for long call or long
% put, reflecting the fact that there are time values in options and as
% time goes, this time value decays. Note particularly that the gap between
% two curves is the most significant (relatively speaking) when S is near
% the strike, this is because at the far left or the far right, the time
% value of the option is almost negligible. On th far right, the value of
% the option is almost 0, and we expect 0 payoff. On the far left, it might
% best to exercise right now. (not sure)
%
% For figure(2) delta plot, we see a more complicated relationship. First
% of all, the delta are all nonpositive. This is because as S increases,
% our expected payoff will decrease, since we are holding power put. Note
% that two curves intersect, let's denote the intersection by x* (around
% 56). Also, we consider absolute value of delta in the following
% comparison of size. So we say size of Delta 1 is larger than that of
% Delta 2, we actually mean that the absolute value of Delta 1 is larger
% than absolute value of Delta 2. Before x*, the size delta for t = 0 is
% less than that of t = T/2. This makes sense since at T/2, any increase in
% asset value at region that already has low asset value (i.e. high payoff)
% will very likely result in a decrease of final payoff, whereas at t =0,
% that phenomenon is less severe. (not sure about this interpretation)
%
% For figure(3) gamma plot, we see a more complicated relationship. First
% of all, the gamma are all nonpositive, indicating that the option value
% is convex and delta is increasing function of S.  Note
% that two curves intersect, let's denote the intersection by x** (around
% 70). Before x**, option value for t = T/2 is more convex. After x**,
% option value for t =0 is more convex.
%
% I am not sure how can I intuitively or financially interpret these plots.

function out = sigmalvf(S)
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
    alpha = 15;
    S0 = 76;
    out = alpha./ (max(S0 + (S- S0)./ (2*S0), 10^(-8)));
end