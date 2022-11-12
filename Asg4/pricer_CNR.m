%% Q6a, only CNR, for different stock, both powerput and standard call and put

function [V1,V2,V3,Vpp] = pricer_CNR(S, N, opt, Tlist, Klist, r)
% output option value at t =0, for all gridpoints, using CN-Rannacher
%
% if opt = 1, compute call option prices and compute powerput on stock A,
% the order of computation matches the order of Tlist, Klist 
% V1, V2, V3 will be option values on the grid S for the 3 call options 
% 
% if opt is sth else, compute put option prices and compute powerput 
% on stock B  
% V1, V2, V3 will be option values on the grid S for the 3 put options 
%
% in both cases, Vpp are the option values of the grid for the powerput

m = length(S);


%% preprocessing step for simple BS equation
% so that we have positive coefficient discretization
% these will also be used for CN and CN-Rannacher
Smid = S(2:end-1);
Sleft = S(1:end-2);
Sright = S(3:end);

if opt == 1
    sigma = sigmalvf_a(Smid);
else
    sigma = sigmalvf_b(Smid);
end

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

%% Boundary condition concern: 
% initialize payoff, call option payoff
% Note that we adopt just the payoff as the boundary condition. By p.9 in
% Lec16, it is suggested to use S_m for call option. But I think the
% difference is small and adopt S_m - K instead, since S_m is quite large.
% similar approach is taken for put and powerput 

%% main for loop 
result= zeros(m,4);

for i = 1:4
    %% Fully implicit
    T = Tlist(i);
    Deltatau = T/N;
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

    %% Crank-Nicolson
    Bin = Deltatau/2* [[[-alphavec;0];0], [[r; (alphavec + betavec +r)];0],[[0;0];-betavec]];
    Mhat = spdiags(Bin, -1:1, m, m);
    % Note that our Mhat is time independent
    I = speye(m);
    [L1,U1,P1] = lu(Mhat+I); 

    %% CN-Rannacher
    K = Klist(i);
    Vold = zeros(m,1); 
    Vnew = Vold;
    % initialize payoff, call option payoff
    if opt == 1
        if i < 4
            Vold = max(S-K,0);
        else 
            Vold = max(K-S,0).^3;
        end
    else 
        if i < 4
            Vold = max(K-S,0);
        else 
            Vold =  max(K-S,0).^3;
        end
    end
        
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
    result(:,i) = Vnew;
end

V1 = result(:,1);
V2 = result(:,2);
V3 = result(:,3);
Vpp =result(:,4);


end
