%% Cleaning the workspace and the command window and adjusting settings
clear;clc
% changeCobraSolver('gurobi', 'LP');

%% Predict usage of underground reations
% model = loadEcModel('../models/ecKmarx.yml');

load('../models/ecKmarx.mat');
% load('../models/Kluyveromyces_marxianus-GEM.mat');

%% Solve batch simulation
modelBatch = ecModel;
% modelBatch = model;

% Carbon sources
% modelBatch = changeRxnBounds(modelBatch, 'r_1726', -13.71, 'b'); % glucose

modelBatch = changeRxnBounds(modelBatch, 'r_1726', 0, 'b'); % glucose

% modelBatch = changeRxnBounds(modelBatch, 'r_1784', -4.71802841918295, 'b'); % fructose
% modelBatch = changeRxnBounds(modelBatch, 'r_1785', -3.66349164057817, 'b'); % galactose
% modelBatch = changeRxnBounds(modelBatch, 'r_1856', -2.19106047326906, 'b'); % lactose
% modelBatch = changeRxnBounds(modelBatch, 'r_1890', -2.19108607754692, 'b'); % sucrose
modelBatch = changeRxnBounds(modelBatch, 'r_1795', -1.38, 'b'); % xylose

modelBatch = changeRxnBounds(modelBatch, 'r_1759', 0, 'b');
modelBatch = changeRxnBounds(modelBatch, 'r_1758', 0, 'b');
modelBatch = changeRxnBounds(modelBatch, 'r_1878', 0, 'b');

modelBatch = changeRxnBounds(modelBatch, 'prot_pool_exchange', -270, 'l');

solutionBatch = optimizeCbModel(modelBatch);

printFluxes(modelBatch, solutionBatch.x, true);
% printFluxes(modelBatch, solutionBatch.x, false);

%% Simulate chemostat (and respiro-fermentative metabolism)
% Code adapted from the plotCrabtree.m function from GECKO3
% model = loadEcModel('../models/ecKmarx.yml');
modelChemo = ecModel;

gRate = 0:0.025:0.475;

results  = zeros(numel(modelChemo.rxns),numel(gRate));

carbon = getIndexes(modelChemo,'r_1726','rxns');
modelChemo = changeObjective(modelChemo, 'r_1913', 0);
modelChemo = changeObjective(modelChemo, 'r_1726', 1);

modelChemo = changeRxnBounds(modelChemo, 'prot_pool_exchange', -200, 'l');
ptot  = -modelChemo.lb(strcmp(modelChemo.rxns,'prot_pool_exchange'));

for i=1:numel(gRate)
    tempModel = modelChemo;

    tempModel = changeRxnBounds(tempModel, 'r_1759', 0, 'b');
    tempModel = changeRxnBounds(tempModel, 'r_1758', 0, 'b');
    tempModel = changeRxnBounds(tempModel, 'r_1878', 0, 'b');

    tempModel = changeRxnBounds(tempModel, 'r_1913', gRate(i), 'l');
    solution = optimizeCbModel(tempModel);
    % printFluxes(tempModel, solution.x, true);

    if ~isempty(solution.x)
        tempModel = changeRxnBounds(tempModel, 'r_1726', solution.x(carbon)*1.05, 'l');
        tempModel = changeRxnBounds(tempModel, 'r_1726', solution.x(carbon)*0.95, 'u');
        tempModel = changeObjective(tempModel, 'r_1726', 0);

        enzymes=find(~cellfun('isempty',strfind(tempModel.rxnNames,'prot_')));
        enzNames=tempModel.rxnNames(enzymes); 
        tempModel=changeObjective(tempModel,enzNames,1);

        solution2 = optimizeCbModel(tempModel);
        % printFluxes(tempModel, solution2.x, true);

        results(:,i) = solution2.x;
    end
end

%%
rxnsToCheck = getIndexes(modelChemo,{'r_1725','r_1775','r_1726','r_1803'},'rxns');

tiledlayout(1,2)

% Plot fluxes
nexttile
plot(gRate,abs(results(rxnsToCheck,:)));
hold on
legend(modelChemo.rxnNames(rxnsToCheck),'Location','northwest')
xlabel('Growth rate (/hour)')
% ylim([0 20])
ylabel('Absolute flux (mmol/gDCWh)')
hold off

% Plot total protein usage
nexttile
poolRxn = getIndexes(modelChemo,'prot_pool_exchange','rxns');
plot(gRate,abs(results(poolRxn,:))/ptot)
xlabel('Growth rate (/hour)')
ylabel('Fraction of protein pool used')
fig=gcf;
fig.Position(3:4) = [700 300];
