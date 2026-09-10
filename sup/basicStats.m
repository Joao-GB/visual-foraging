%% One-way ANOVA
% y = response
% group = categorical group indicator

[p,tbl,stats] = anova1(y, group);

% Post-hoc comparisons
multcompare(stats);

%% Linear regression
% x = continuous predictor
% y = response

mdl = fitlm(x, y);

disp(mdl)
plot(mdl)

%% Continous + categorical predictor
% x = continuous predictor
% group = categorical
% y = response

T = table(x, group, y);

mdl = fitlm(T, 'y ~ x + group');

disp(mdl)

% Comparar duas regressões (?)
T = table(x, group, y);

mdl = fitlm(T, 'y ~ x * group');

disp(mdl)
plotInteraction(mdl, 'x', 'group')

%% Relação não linear
T.x2 = T.x.^2;

mdl = fitlm(T, 'y ~ x + x2 + group + x:group + x2:group');

disp(mdl)