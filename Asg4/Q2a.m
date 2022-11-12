%% Q2a, variable timestepping and constant timestepping

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

% make convergence table
CNRvartable = zeros(length(Deltataulist),7);
CNRconstanttable = zeros(length(Deltataulist),7);
vartimetotal =0;
contimetotal =0;
for i = 1:5
    % obtain option value at t=0, S = S0 using the corresponding grid and
    % time step
    Deltatau = Deltataulist(i);
    S = Scell{i};
    dnorm = dnormlist(i);
    
    
    vartime = tic;
    [V0var, steptotalvar, countervar] = CNR(S, Deltatau, S0, T, K, r,dnorm,1);
    vartime = toc(vartime);
    vartimetotal = vartimetotal + vartime;
    
    contime = tic;
    [V0con, steptotalcon, countercon] = CNR(S, Deltatau, S0, T, K, r,dnorm,0);
    contime = toc(contime);
    contimetotal = contimetotal + contime;

    % nodes number 
    CNRvartable(i,1) = size(Scell{i},1);
    CNRconstanttable(i,1) = size(Scell{i},1);

    % timesteps number
    CNRvartable(i,2) = T/ Deltatau;
    CNRconstanttable(i,2) = T/ Deltatau;

    % option Value
    CNRvartable(i,3) = V0var;
    CNRconstanttable(i,3) = V0con;
    % change in option value 
    if i > 1
        CNRvartable(i,4) = CNRvartable(i,3) - CNRvartable(i-1,3);
        CNRconstanttable(i,4) = CNRconstanttable(i,3) - CNRconstanttable(i-1,3);
    end
    
    % Ratio
    if i > 2 
        CNRvartable(i,5) = CNRvartable(i-1,4)/CNRvartable(i,4);
        CNRconstanttable(i,5) = CNRconstanttable(i-1,4)/CNRconstanttable(i,4);
    end
    
    % counter 
    CNRvartable(i,6) = countervar;
    CNRconstanttable(i,6) = countercon;

    % total iterative steps
    CNRvartable(i,7) = steptotalvar;
    CNRconstanttable(i,7) = steptotalcon;
    
end

CNRvartable = array2table(CNRvartable, 'VariableNames', ...
    {'Nodes','initial timestep','Value','Change','Ratio','# timestepping','iterative steps' })
vartimetotal

CNRconstanttable = array2table(CNRconstanttable, 'VariableNames', ...
    {'Nodes','timesteps','Value','Change','Ratio','# timestepping','iterative steps' })
contimetotal

%% Discussion
% For both methods, we refine the grids as in Q1, doubling the number of
% intervals.
%
% For CN-Rannacher with variable timestepping, at each refinement, we
% reduce the initial timestep by a factor of 4, and reduce dnorm by a
% factor of 2. 
% 
% For CN-Rannacher with constant timesteppin, at each refinement, we reduce
% the constant timestepping by a factor of 4. 
% 
% We see that CN-Rannacher with variable timestepping does give us
% quardartic convergence, since the ratio in CNRvartable is around 4. 
% This is expected.
% 
% Note that for the convergence table of CN-Rannacher with constant
% timestepping, at each refinement we are dividing the constant timestep by
% 4. This might explain why we also observe quardartic convergence in such
% case. I suspect that the (convergence) error might only depend slightly
% on \Delta S, so that the number of nodes does not greatly affect the
% convergence rate. If this is true, this explains why we obtain quardartic
% convergence. Another reason for that might be the powerput payoff is
% quite smooth so that quardartic convergence can still be achieved.
% 
% Note that if we change the refinement to timestep /2, still observe the
% quardartic convergence... (WHY???) Perhaps it is because the powerput
% payoff is quite smooth. I have experimented my codes on other examples in
% the lecture slides and obtain similar results.
%
% Lastly, in the setting we consider here, the runtime of CNR-variable
% timestepping is less than that of CNR-constant timestepping.

function [V0, steptotal, counter] = CNR(S, Deltatau, S0, T, K, r,dnorm,opt)
%% Description 
% Inputs:
% opt = 1 for variable time, else for constant timestepping
% when opt = 1, Deltatau is the initial timestep 
% when opt is else, Deltatau is the constant timestep 
% Outputs:
% V0, option value at S0 at t = 0 
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
        dtnew = (dnorm/maxrelchange) * dtold;
        dtold = dtnew;
        % update V
        Vold = Vnew;
    end

    S0id = S0 == S;
    V0 = Vnew(S0id);
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
    S0id = S0 == S;
    V0 = Vnew(S0id);
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

% why
% largediagonal = zeros(m,1);
% largeidx = Vnew < Vstar;
% largediagonal(largeidx) = Large;
% P = spdiags(largediagonal, 0,m,m);
% sum(P,'all')
end

function out = sigmalvf(S)
% note that the local volatility function actually is independent of t 
% vectorized output, assume column vector input 
    alpha = 15;
    S0 = 76;
    out = alpha./ (max(S0 + (S- S0)./ (2*S0), 10^(-8)));
end