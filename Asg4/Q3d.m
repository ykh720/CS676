%% Q3d

clear vars
close all

%% Data preparation 
S0 =100;
r = 0.03;

Klist = [0.9, 0.95, 1, 1.05, 1.1]* S0;
Tlist = [0.425, 0.695, 1];

% IV table
table2 = zeros(3,5);
table2(1,:) = [0.155, 0.138, 0.125, 0.109, 0.103];
table2(2,:) = [0.157, 0.144, 0.133, 0.118, 0.104];
table2(3,:) = [0.159, 0.149, 0.137, 0.127, 0.113];

market = zeros(3,5);

for j = 1:5
    for i = 1:3
        [Call , put] = blsprice(S0, Klist(j), r, Tlist(i), table2(i,j));
        market(i,j) = Call;
    end
end

market = reshape(market, numel(market), []);

%% Levenberg-Marquardt method
options = optimset('Jacobian', 'on', 'Algorithm', 'levenberg-marquardt',...
    'Display', 'iter', 'MaxIter', 50);

% initial x 
x0 = [0,0,1,1,0,0]';

% xopt: optimal parameters xopt
% cal_error: calibration error
[xopt, cal_error] = lsqnonlin(@myfun, x0, [],[], options) 

% for more information on the iterative display, see 
% https://www.mathworks.com/help/optim/ug/iterative-display.html#f92387

%% IV from true model and the IV from the optimal parameters 
Vfromxopt = myfun(xopt) + market;
Vfromxopt = reshape(Vfromxopt, length(Tlist), length(Klist))

[X,Y] = meshgrid(Klist, Tlist);
IVfromxopt = blsimpv(S0, X, r, Y, Vfromxopt)

% IVtrue = blsimpv(S0, X,r,Y, reshape(market, 3,5))
% note that table 2 are the implied volatilties from the true model
table2

figure(1)
mesh(X,Y, table2,'FaceColor','blue')
xlabel('K')
ylabel('T')
zlabel('IV')
title('implied volatility surface')
hold on 
mesh(X,Y, IVfromxopt,'FaceColor','red')
legend('true','from xpot')

relerrorIV = max(abs(IVfromxopt- table2)./table2,[],'all')

%% Discussion
% We see that in general, our implied volatilties calibration does a fairly
% good job, with the relative error being around 12%. That is, the implied
% volatilties from the option values computed by model as given by xopt,
% are roughly 88% accurate, when compared to the true model, i.e. the
% implied volatilities in table 2. Based on the simplicity of our model
% (simple neural network), this is already quite good.
% 
% Some thoughts on the formulation of the minimization problem. Note that
% the objective function is the 1/2 * sum of square of difference between
% option values from our model and option values from true model (or
% so-called market prices). What if we directly incorporate the difference
% between implied volatilties in the objective function? For instance, we
% can define F to be the sum of the square of the difference between the
% true implied volatilities and the implied volatilities from our model.
% What are the advantages and disadvantages of such method?