%% Q3b, only CNR 

function V0CNRan = callpricer_CNR(S, Deltatau, x, T, K, r)
% output option value at t =0, for all gridpoints, using fully implicit,
% Crank Nicolson, and CN-Rannacher
%
% input: x is the set of parameters. x =(w_1,w_2,u_1,u_2,v_1,v_2)', used in
% local volatility

m = length(S);
Vold = zeros(m,1); 
Vnew = Vold;
N = T / Deltatau;
%% Boundary condition concern: 
% initialize payoff, call option payoff
% Note that we adopt just the payoff as the boundary condition. By p.9 in
% Lec16, it is suggested to use S_m for call option. But I think the
% difference is small and adopt S_m - K instead, since S_m is quite large.
Vold = max(S-K,0);

%% preprocessing step for simple BS equation
% so that we have positive coefficient discretization
% these will also be used for CN and CN-Rannacher
Smid = S(2:end-1);
Sleft = S(1:end-2);
Sright = S(3:end);
sigma = sigmalvf(Smid,x);
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

%% Fully implicit

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

% for n = 1 : N  
%     y = L\(P*Vold);
%     Vnew = U\y;
%     Vold = Vnew;
% end
% V0fullyimplicit = Vnew;


%% Crank-Nicolson
%%% only compute initial option value
%%% can easily modified to obtain options values at each grid points
% Vold = zeros(m,1); 
% Vnew = Vold;
% % initialize payoff, call option payoff
% Vold = max(S-K,0);

Bin = Deltatau/2* [[[-alphavec;0];0], [[r; (alphavec + betavec +r)];0],[[0;0];-betavec]];
Mhat = spdiags(Bin, -1:1, m, m);
% Note that our Mhat is time independent
I = speye(m);
[L1,U1,P1] = lu(Mhat+I); 


% for n = 1 : N  
%     y = L1\(P1*(I-Mhat)*Vold);
%     Vnew = U1\y;
%     Vold = Vnew;
% ends
% V0CN = Vnew;

%% CN-Rannacher
Vold = zeros(m,1); 
Vnew = Vold;
% initialize payoff, call option payoff
Vold = max(S-K,0);

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
V0CNRan = Vnew;

end

function out = sigmalvf(S,x)
%%  local volatility by simple neural network, (9) in asg4
% 
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
%
% inputs: S, gridpoints of S 
%         x, set of parameters. x =(w_1,w_2,u_1,u_2,v_1,v_2)'
w1 = x(1);
w2 = x(2);
u1 = x(3);
u2 = x(4);
v1 = x(5);
v2 = x(6);
y1 = 1./(u1 + v1 * S);
y2 = 1./(u2 + v2 *S);
out = 1./(1 + exp(w1 * y1 + w2* y2));
end