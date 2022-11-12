function [V0, steptotal, counter] = CNRall(S, Deltatau, S0, T, K, r,dnorm,opt)
%% Description 
% Inputs:
% opt = 1 for variable time, else for constant timestepping
% when opt = 1, Deltatau is the initial timestep 
% when opt is else, Deltatau is the constant timestep 
% Outputs:
% V0, option value at S (all the grid points) at t = 0 
% steptotal, the total number of steps used for iterative solvers
% counter, the total time of choosing variable timestep; counter is just
% T/Deltatau when constant timestepping 

%%  preprocessing step for simple BS equation
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
%% Some common data needed for constant and variable timestepping 
m = length(S);
Vold = zeros(m,1); 
Vnew = Vold;
N = T / Deltatau;
% initialize payoff 
Vold = max(K-S,0).^3;

% penalty method data, tolerance and Large (penalty coefficient)
tol = 10^(-6); % as suggested by piazza 
Large = 1/tol;

% Vstar, payoff at the gridpoints, for American option payoff constraints
Vstar = max(K-S,0).^3;

Bin = [[[-alphavec;0];0], [[r; (alphavec + betavec +r)];0],[[0;0];-betavec]]; 
% generic bin, used repeatedly below 

steptotal = 0;

%% opt = 1, Use variable timestep 
%%% also obtain number of steps to reach tol
%%% also need to update the Mhat and M 
if opt == 1     
    dtold = Deltatau;
    dtnew = 0;
    currenttau = 0;
    D = 1; % p.32 in lec 18, D = 1 fine for dollars
    counter = 0;
    % CN-Rannacher, variable time
    while currenttau < T
        counter = counter +1;

        if currenttau + dtold >= T
            dtold = T - currenttau;
            currenttau = T;
        else
            currenttau = currenttau + dtold;
        end

        if counter <= 2 % first two steps are fully implicit
            % build M, using correct \Delta tau
            M = spdiags(dtold*Bin, -1:1, m, m);
            [Vnew,step]= iterativesolver(Vold, Vstar, tol, Large, M,m,0);
        else 
            % build Mhat, using correct \Delta tau
            Mhat = spdiags((dtold/2) *Bin, -1:1, m, m);
            [Vnew,step]= iterativesolver(Vold, Vstar, tol, Large, Mhat,m,1);
        end

        steptotal = steptotal + step;

        % chooose timestep at the end
        % assume Vold is still storing previous data, not current data
        inter = max(D,max(abs(Vnew), abs(Vold)));
        maxrelchange = max(abs(Vnew - Vold)./inter);
        dtnew = dnorm/maxrelchange * dtold;
        dtold = dtnew;
        % update V
        Vold = Vnew;
    end

    V0 = Vnew;
else
    %% opt is not 1, use constant timestepping
    Mhat = spdiags((Deltatau/2)*Bin, -1:1, m, m);
    M = spdiags(Deltatau*Bin, -1:1, m, m);
    % two steps of Fully Implicit
    for i = 1:2
        [Vnew,step]= iterativesolver(Vold, Vstar, tol, Large, M,m,0);
        steptotal = steptotal + step;
        Vold = Vnew;
    end
    
    for i = 1: N-2
        [Vnew,step]= iterativesolver(Vold, Vstar, tol, Large, Mhat,m,1);
        steptotal = steptotal + step;
        Vold = Vnew;
    end
    V0 = Vnew;
    counter = N;
end

end

function [Vnew,step]= iterativesolver(Vold, Vstar, tol, Large, genM,m,opt)
% opt = 1 for CN, else fully implicit
% genM= Mhat when opt = 1 (i.e. CN), o/w genM = M (i.e. fully implicit)

error = 10000; % set the error to be some large number 
step = 0;
Viterold = Vold;
I = speye(m);
while error >= tol
    % Vold is V^n
    % Viterold is (V^{n+1})^k, Viternew will be (V^{n+1})^(k+1)
    largediagonal = zeros(m,1);
    largeidx = Viterold < Vstar;
    largediagonal(largeidx) = Large;
    P = spdiags(largediagonal, 0,m,m);
    if opt == 1 % CN
        y = (I - genM)*Vold + P* Vstar;
        Viternew = (I + genM + P) \ y;
    else % fully implicit
        y = Vold + P* Vstar;
        Viternew = (I + genM + P) \ y;
    end
    
    inter = abs(Viternew - Viterold)./max(1,abs(Viternew));
    error = max(inter); 
    step = step + 1;
    if error < tol
        Vnew = Viternew;
        break;
    end
    Viterold = Viternew;
end

end

function out = sigmalvf(S)
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
    alpha = 15;
    S0 = 76;
    out = alpha./ (max(S0 + (S- S0)./ (2*S0), 10^(-8)));
end