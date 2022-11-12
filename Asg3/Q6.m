%% Q6a
clearvars
close all

load('RawData.mat')

% change them all to column vectors
CSTest = CSTest';
CSTrain = CSTrain';
CVTest = CVTest';
CVTrain = CVTrain';

%% hedging error for training set 
% assume the data is sorted with ascending time 

dailyhedgingerror_training = CVTrain - DeltaTrain.* CSTrain;

g(1) = figure(1);
histogram(dailyhedgingerror_training,50)
title('histogram for daily hedging error, BS, training')

resulttable= zeros(1,4);
resulttable(1) = mean(dailyhedgingerror_training);
resulttable(2) = std(dailyhedgingerror_training);
beta = 0.95;
[var,cvar] = dVaRCVaR(dailyhedgingerror_training, beta);
resulttable(3) = var;
resulttable(4) = cvar;
resulttable_training_BS = array2table(resulttable, 'VariableNames',{'mean',...
    'standard deviation', '95% VaR','95% CVaR'});

%% hedging error for testing set 
% assume the data is sorted with ascending time 

dailyhedgingerror_testing = CVTest - DeltaTest.* CSTest;

g(2) = figure(2);
histogram(dailyhedgingerror_testing,50)
title('histogram for daily hedging error, BS, testing')

resulttable= zeros(1,4);
resulttable(1) = mean(dailyhedgingerror_testing);
resulttable(2) = std(dailyhedgingerror_testing);
beta = 0.95;
[var,cvar] = dVaRCVaR(dailyhedgingerror_testing, beta);
resulttable(3) = var;
resulttable(4) = cvar;
resulttable_testing_BS = array2table(resulttable, 'VariableNames',{'mean',...
    'standard deviation', '95% VaR','95% CVaR'});

%% Q6b 1
% we use OLS linear regression tot learn the parameters a,b,c 
% Note that we use (4) \Delta f = \delta_{MV} \Delta S + \epsilon, in the
% paper and make use of (7) in the assignment. 
% After some algebraic manipulation, we obtain the following form for
% linear regression
% \Delta f - \delta_{BS} \Delta S = \Delta S Vega/(S\sqrt{T}) a + 
% \Delta S Vega \delta_{BS} /(S\sqrt{T}) b 
% + \Delta S Vega \delta_{BS}^2 /(S\sqrt{T}) c + \epsilon
%
% Note that \Delta f - \delta_{BS}\Delta S is the BS hedging error
%
% Use the training data
inter = CSTrain.* VegaTrain./(STrain.* sqrt(TauTrain));
X = [ inter, inter.* DeltaTrain, inter.*DeltaTrain.*DeltaTrain];
mdl = fitlm(X, dailyhedgingerror_training,'Intercept', false)
a = mdl.Coefficients{1,1};
b = mdl.Coefficients{2,1};
c = mdl.Coefficients{3,1};

delta_MV_training = DeltaTrain + (VegaTrain./(STrain.* sqrt(TauTrain))).*...
    (a + b * DeltaTrain + c * DeltaTrain.^2);

% alternatively, can use anova function to get sum of squared error
anova(mdl,'summary');

SSE_MV_training = sum((CVTrain- CSTrain .* delta_MV_training).^2);

Gain_training = 1 - SSE_MV_training ./ ( sum(dailyhedgingerror_training.^2));

delta_MV_testing = DeltaTest + (VegaTest./(STest.* sqrt(TauTest))).*...
    (a + b * DeltaTest + c * DeltaTest.^2);

SSE_MV_testing = sum((CVTest- CSTest .* delta_MV_testing).^2);

Gain_testing = 1 - SSE_MV_testing ./ ( sum(dailyhedgingerror_testing.^2));

%% Q6b 2
hedgeerrorMV_training = CVTrain- CSTrain .* delta_MV_training;
g(3) = figure(3);
histogram(hedgeerrorMV_training,50)
title('histogram for daily hedging error, MV, training')

resulttable= zeros(1,4);
resulttable(1) = mean(hedgeerrorMV_training);
resulttable(2) = std(hedgeerrorMV_training);
beta = 0.95;
[var,cvar] = dVaRCVaR(hedgeerrorMV_training, beta);
resulttable(3) = var;
resulttable(4) = cvar;
resulttable_training_MV = array2table(resulttable, 'VariableNames',{'mean',...
    'standard deviation', '95% VaR','95% CVaR'});


hedgeerrorMV_testing = CVTest- CSTest .* delta_MV_testing;
g(4) = figure(4);
histogram(hedgeerrorMV_testing,50)
title('histogram for daily hedging error, MV, testing')

resulttable= zeros(1,4);
resulttable(1) = mean(hedgeerrorMV_testing);
resulttable(2) = std(hedgeerrorMV_testing);
beta = 0.95;
[var,cvar] = dVaRCVaR(hedgeerrorMV_testing, beta);
resulttable(3) = var;
resulttable(4) = cvar;
resulttable_testing_MV = array2table(resulttable, 'VariableNames',{'mean',...
    'standard deviation', '95% VaR','95% CVaR'});

% for easy comparsion
resulttable_training_BS;
resulttable_testing_BS;

%% Discussion Q6b 2 
% We compare MV delta on training set against BS delta on training set and
% compare MV delta on testing set against BS delta on testing set. From the
% histograms, we see that the histograms for MV delta are more concentrated
% around 0 since the x-axis range of it is smaller. Moreover, by the table
% (mean, standard deviation, VaR, CVaR), we see that MV delta does give
% smaller standard deviation. This can be seen by comparing
% resulttable_training_BS vs resulttable_training_MV and comparing
% resulttable_testing_BS vs resulttable_testing_MV. We also see that by
% employing MV, we have improvments in VaR and CVaR, both in training and
% testing sets, compared to using BS delta.
%
% Now we compare MV delta on training set against MV delta on testing set.
% We see that the mean on training set is larger than that on testing set.
% Such phenomenon also exists on BS delta case. But in our case the
% difference between the mean is actually lower than that in BS delta case.
% More importantly, we observe that the standard deviation, VaR, CVaR (in
% absolute value) for training set is lower than that of testing set. That
% means hedging performance on training set is better. This makes sense
% since the coefficients a,b,c are fitted using training set so we expect
% better behaviour on training set. Keep in mind that we do have improvment
% on testing case, as compared to using BS delta, as reflected by
% Gain_testing being around 0.34.

%% Q6c

% design matrix
Xnew = [ones(length(STrain),1), STrain, DeltaTrain, DeltaTrain.^2, VegaTrain .* DeltaTrain, ...,
    VegaTrain.^2, DeltaTrain.^3].* CSTrain;
mdl_c = fitlm(Xnew, CVTrain, 'Intercept',false)

test_coeff=  regress(CVTrain, Xnew);

% Note that the p-values of c4 and c6 are pretty big, that means they might
% not be significant variables.

c0 = mdl_c.Coefficients{1,1};
c1 = mdl_c.Coefficients{2,1};
c2 = mdl_c.Coefficients{3,1};
c3 = mdl_c.Coefficients{4,1};
c4 = mdl_c.Coefficients{5,1};
c5 = mdl_c.Coefficients{6,1};
c6 = mdl_c.Coefficients{7,1};

% b 1

delta_c_train = c0 + c1* STrain + c2*DeltaTrain + c3 * DeltaTrain.^2 + ...
    c4 * VegaTrain .* DeltaTrain + c5 * VegaTrain.^2 + c6* DeltaTrain.^3;

SSV_c_train = sum((CVTrain- CSTrain .* delta_c_train).^2);

Gain_c_train = 1 - SSV_c_train ./ ( sum(dailyhedgingerror_training.^2));

delta_c_test = c0 + c1* STest + c2*DeltaTest + c3 * DeltaTest.^2 + ...
    c4 * VegaTest .* DeltaTest + c5 * VegaTest.^2 + c6* DeltaTest.^3;

SSE_c_testing = sum((CVTest- CSTest .* delta_c_test).^2);

Gain_c_test = 1 - SSE_c_testing ./ ( sum(dailyhedgingerror_testing.^2));

% b 2 

hedgeerror_c_train = CVTrain- CSTrain .* delta_c_train;
g(5) = figure(5);
histogram(hedgeerror_c_train,50)
title('histogram for daily hedging error, part c, training')

resulttable= zeros(1,4);
resulttable(1) = mean(hedgeerror_c_train);
resulttable(2) = std(hedgeerror_c_train);
beta = 0.95;
[var,cvar] = dVaRCVaR(hedgeerror_c_train, beta);
resulttable(3) = var;
resulttable(4) = cvar;
resulttable_train_c = array2table(resulttable, 'VariableNames',{'mean',...
    'standard deviation', '95% VaR','95% CVaR'});

hedgeerror_c_test = CVTest- CSTest .* delta_c_test;
g(6) = figure(6);
histogram(hedgeerror_c_test,50)
title('histogram for daily hedging error, part c, testing')

resulttable= zeros(1,4);
resulttable(1) = mean(hedgeerror_c_test);
resulttable(2) = std(hedgeerror_c_test);
beta = 0.95;
[var,cvar] = dVaRCVaR(hedgeerror_c_test, beta);
resulttable(3) = var;
resulttable(4) = cvar;
resulttable_test_c = array2table(resulttable, 'VariableNames',{'mean',...
    'standard deviation', '95% VaR','95% CVaR'});

%% Result
% we present all the results here cleanly
result = [resulttable_training_BS;...
    resulttable_testing_BS; resulttable_training_MV;...
    resulttable_testing_MV;resulttable_train_c; resulttable_test_c];
result.Properties.RowNames = {'BS,train','BS,test','MV,train','MS,test',...
    'part c,train','part c,test'};
result 

Gain = zeros(1,4);
Gain(1) = Gain_training;
Gain(2) = Gain_testing;
Gain(3) = Gain_c_train;
Gain(4) = Gain_c_test;

Gain = array2table(Gain, 'VariableNames',{'MV,train','MV,test','part c,train',...
    'part c,test'})

for i =1:length(g)
    saveas(g(i),strcat('fig_Q6_',string(i)),'epsc')
end

%% Q6c Discussion 
% Note that generally part c parametric form behaves quite well and achieve
% similar improvements as in MV delta. It is because the parametric form
% here is quite powerful, we included a lot of predictors and allow up to
% cubic delta. In a regression point of view, more predictors will give
% better training loss, while susceptible to worse generalization loss,
% especically in the case of adding irrelevant predictors. Observing the
% p-values of the predictors in part c, the p-value of square of vega is
% quite large meaning that it might not be relevant. In fact, the added
% complexity in this part c parametric form results into behaving better on
% training set while behaving worse on testing set, as compared to MV
% delta. This can be seen by observing the standard deviation and the Gain.
% Note that Gain on testing for part c is lower than that from MV delta.
% Note that Gain on training for part c is higher than that from MV delta.
% This observation is reasonable since our complex model will fit well on
% training set (overfitting) but not necessarily fit well on testing set.
% One reason for delta MV behaves better on the testing set is that it is
% supported by empirical observations of S&P 500 options on the
% relationship of delta_MV and delta_BS. Note in particular that there is
% no inverse of (S sqrt (T)) term in part c. That lack of domain knowledge
% might explain the worse behavior of part c on testing set and in
% generalization.
%
% But after all, generally the parametric form of part c captures quite a
% lot information on MV delta so that the resulting standard deviation and
% sum of square error than just using BS delta.