%% Q6a 
clearvars
close all

rng('default')

SA0scalar = 20;
SB0scalar = 25;
r = 0.03;
m = 10000;

SA0 = ones(m,1) * SA0scalar;
SB0 = ones(m,1) * SB0scalar;

T = 1/12; % one month 
N = 10;

SA = eulerlvf(SA0, T, N, 1, r);
SA = sort(SA); % sort the SA so that the order matches options 
% max(SA) -20 
% min(SA)- 20
SB = eulerlvf(SB0, T, N, 0, r);
SB = sort(SB);
% max(SB) - 25
% min(SB) -25 

SAreturn = (SA - SA0scalar) ./ SA0scalar;
SBreturn = (SB - SB0scalar) ./ SB0scalar;

% make grid, refine the grid 4 times 
SAgrid = gridrefine(SA0scalar, SA,4); 
SBgrid = gridrefine(SB0scalar, SB,4);

% note that the probability of obtaining duplicate elements are 0 
% in fact our choice of random seed make sure that there is no duplicate in
% SAgrid, SBgrid

% save all the rate of return on all the instruments 
return_A = zeros(m,5);
return_B = zeros(m,5);
% first column stores the stock
return_A(:,1) = SAreturn;
return_B(:,1) = SBreturn;

Tlist = [0.5,0.5,0.5, 2/3]';
Klist = [SA0scalar, 1.1*SA0scalar, 1.2*SA0scalar, SA0scalar]';
[~,idxA] = intersect(SAgrid, SA);
N = 100;
idscalar = SAgrid == SA0scalar;
[V1old,V2old,V3old,Vppold] = pricer_CNR(SAgrid, N, 1, Tlist, Klist, r);
price_old = [V1old(idscalar),V2old(idscalar),V3old(idscalar),Vppold(idscalar)];
[V1,V2,V3,Vpp] = pricer_CNR(SAgrid, N, 1, Tlist- 1/12, Klist, r);
price_new = [V1(idxA),V2(idxA),V3(idxA),Vpp(idxA)];
return_A(:,2:5) = (price_new - price_old) ./ price_old;
% for stock A, 2nd to 4th columns are call, last column is powerput 

Tlist = [1,1,1,1]';
Klist = [0.8*SB0scalar, 0.9*SB0scalar, SA0scalar, SA0scalar]';
[~,idxB] = intersect(SBgrid, SB);
N = 100;
idscalar = SBgrid == SB0scalar;
[V1old,V2old,V3old,Vppold] = pricer_CNR(SBgrid, N, 0, Tlist, Klist, r);
price_old = [V1old(idscalar),V2old(idscalar),V3old(idscalar),Vppold(idscalar)];
[V1,V2,V3,Vpp] = pricer_CNR(SBgrid, N, 0, Tlist- 1/12, Klist, r);
price_new = [V1(idxB),V2(idxB),V3(idxB),Vpp(idxB)];
return_B(:,2:5) = (price_new - price_old) ./ price_old;
% for stock B, 2nd to 4th columns are put, last column is powerput 


% obtain equally weighted portfolio
eq_weg_portfolio = 1/10 * (sum(return_A, 2) + sum(return_B,2));
% histogram
figure(1)
histogram(eq_weg_portfolio,50)
title('histogram of equal allocation')
xlabel('rate of return')
ylabel('frequency')

% Note that the histogram is not normal, very skewed 

%% Q6b 
allresult = [return_A, return_B]; % 10000 by 10 matrix 

loss = -allresult; % so that we consider a min opt problem  

betalist = [0.8, 0.85, 0.9,0.95];
rholist = zeros(4,1);
M = 10000;
n =10;

avgloss = mean(loss)'; 
f = [avgloss;zeros(M+1,1)]; % so that we consider a min opt problem  
varlength= n + M + 1; % our variable is listed as [x,y,alpha]'
Aeq = zeros(1, varlength);
Aeq(1,1:n) =1;
beq = 1;

A = zeros(n + 1 + M + M,varlength);
b = zeros(n + 1 + M + M,1);
% first n columns for x_i \geq 0 
A(1:n, 1:n) = -eye(n); % need to add negative sign since A x \leq b 
b(1:n) = 0;

% next column for CVar constraint, we will update it in the for loop below 

% then constraints on the variables y that is larger L^T_i x - alpha 
A(n+2:n+1+M,n+1:n+M) = -eye(M);
A(n+2:n+1+M, end) = -1;
A(n+2:n+1+M, 1:n) = loss;
b(n+2:n+1+M) = 0;

% then nonnegative constraint on y 
A(n+1+M+1:end, n+1:n+M) = -eye(M);
b(n+1+M+1:end) = 0;

allvaroptlist = zeros(varlength, 4);

comparetable = zeros(5,5);
alphalist = zeros(4,1 );


for i = 1:4
    beta =betalist(i); % for solving min opt problem
    [~,cvar] = dVaRCVaR(eq_weg_portfolio, beta);
    rholist(i) = -cvar;  % consider cvar in terms of loss
    A(n+1, end) = 1;
    A(n+1,n+1:n+M) = 1/(M*(1-beta));
    b(n+1) = rholist(i);
    % solve min problem
    allvar = linprog(f,A,b,Aeq,beq);
    allvaroptlist(:,i) = allvar;
    x = allvar(1:10);
    alpha = allvar(end);
    
    alphalist(i) = alpha; % var approximation 
    
    optportfolio = allresult * x;
    
%     figure(i+2)
%     histogram(optportfolio,50)
    
    optmean = mean(optportfolio);
    
    comparetable(i+1,1) = optmean;
    
    for j = 1:4 % for table,  
        beta = betalist(j);
        [~,cvar] = dVaRCVaR(optportfolio, beta);
        comparetable(i+1,j+1) = -cvar; % consider cvar in terms of loss
    end
    % note that we can also kind of recover the CVaR by using 
    % alpha + 1 / (M*(1-beta)) * sum(y), but we can only recover one CVaR
    % for the particular beta we used in minimization problem 
end

comparetable(1,1) = mean(eq_weg_portfolio);
comparetable(1,2:end) = rholist;

format short 

comparetable = array2table(comparetable,'VariableNames',...
    {'mean','80% CVaR','85% CVaR','90% CVaR','95% CVaR'},...
    'RowNames',{'equal allocation','optimal wrt 80%', ...
    'optimal wrt 85%','optimal wrt 90%','optimal wrt 95%'})

%% Discussion
% We see that the optimal portfolios (with respect to different beta), have
% very similar mean rate of return, around 0.00363. Note that the rate of
% return we are considering here is for holding the portfolio for one
% month. So the risk-free rate that we should compare to is exp(r/12)-1,
% which is around 0.002503127605795. We observe that the equal allocation
% portfolio has mean rate return 0.00269571202274778, which is very close
% to the risk-free rate. We observe that the mean rate return for optimal
% portfolios are siginificantly better than equal allocation and risk free
% rate.
%
% For optimal portfolio corresponding to beta and rho_beta (i.e CVaR of
% equal allocation portfolio with confidence level beta), the CVaR with
% confidence level beta is very close to rho_beta. This is true because the
% optimization procedure will kind of try to make use of all the given
% budget (i.e. rho_beta) of the CVaR of the optimal portfolio. Note that
% the second line constraint is an upper bound of the CVaR_beta of the
% optimal portfolio. To maximize the return, the optimal portfolio will try
% its best to reach the upper bound, making more rooms for increasing the
% return and hence resulting in a CVaR_beta that is close to rho.
%
% Suppose beta_1 < beta_2. (for instance beta_1 = 0.8, beta_2 = 0.9). We
% observe that optimal portfolio wrt beta_1 has mean slight lower than that
% wrt beta_2. Moreover, with respect to any beta that we consider, the
% CVaR_beta of the former portfolio are always less than the latter
% portfolio. I suspect that it is true because we actually put more
% constraints on the former portfolio since a lower beta_1 means that we
% have to consider more scenarios when calculating the CVaR. Thus, we put
% constraints on more scenarios and resulting in a slightly "safer"
% portfolio than the latter portfolio in the sense that all CVaR calculated
% are smaller. But it is done at the cost of having a lower mean rate of
% return. (in short, constraint on 80% CVaR impose constraints on worse 20%
% scenarios, while constraint on 95% CVaR only impose constraints on worse
% 5% scenarios. A limit on 80% CVaR will often also limit 95% CVaR to a
% satisfiable degree while the converse is not necessarily the case.
%
% Another possible reason for the above phenomenon might be that the risk
% budget obtained from taking CVaR with respect to different beta on the
% equal allocation portfolio acutally result in a rather uneven
% distribution of risk budget that allows more risk tolearance when we take
% large beta. (i think the last paragraph is more convincing)

function S = gridrefine(S0, Snew,n) 
% S0, initial stock value 
% Snew, all new stock prices that we want to add
% n, number of refinements on (2) in asg 4 
S = [ 0:0.1*S0:0.4*S0, ...
    0.45*S0:0.05*S0:0.8*S0, 0.82*S0 :0.02*S0 :0.9 *S0,...
    0.91*S0 :0.01*S0 :1.1*S0, 1.12*S0 :0.02*S0 :1.2*S0, ...
    1.25*S0 :0.05*S0 :1.6*S0, 1.7*S0 :0.1*S0 :2*S0,...
    2.2*S0, 2.4*S0, 2.8*S0, 3.6*S0, 5*S0, 7.5*S0, 10*S0];
S = S';

for i = 1: n
    S = [S;(S(1:end-1) + S(2:end)) / 2 ];
    S = sort(S);
end

S = [S; Snew];
S = sort(S);

end
