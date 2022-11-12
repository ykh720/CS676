%% Q3a 

clear vars
close all

S0 =100;
r = 0.03;

Klist = [0.9, 0.95, 1, 1.05, 1.1]* S0;
Tlist = [0.425, 0.695, 1];

% IV table
table2 = zeros(3,5);
table2(1,:) = [0.155, 0.138, 0.125, 0.109, 0.103];
table2(2,:) = [0.157, 0.144, 0.133, 0.118, 0.104];
table2(3,:) = [0.159, 0.149, 0.137, 0.127, 0.113];

[X,Y] = meshgrid(Klist, Tlist);

figure(1) 
mesh(X,Y, table2)
xlabel('K')
ylabel('T')
zlabel('IV')
title('implied volatility surface')

optiontable = zeros(3,5);

for j = 1:5
    for i = 1:3
        [Call , put] = blsprice(S0, Klist(j), r, Tlist(i), table2(i,j));
        optiontable(i,j) = Call;
    end
end

figure(2) 
mesh(X,Y, optiontable)
xlabel('K')
ylabel('T')
zlabel('option values')
title('option value surface, European Call')

%% Q3a Discussion 
% We see that implied volatility surface and option value surface have
% similar shape. Maybe option value surface is a little bit smoother? 