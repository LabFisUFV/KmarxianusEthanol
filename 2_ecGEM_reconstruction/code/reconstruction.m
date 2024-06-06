%% STAGE 1: expansion from a starting metabolic model to an ecModel structure
% Load adapter and model
adapterLocation = fullfile('../KmarxAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

model = importModel('../models/Kluyveromyces_marxianus-GEM.xml');

% Make ecModel
[ecModel, noUniprot] = makeEcModel(model,false);

saveEcModel(ecModel,'ecKmarx_stage1.yml');

%% STAGE 2: integration of kcat into the ecModel structure
% Gather EC numbers, required for BRENDA as kcat source
ecModel         = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);
ecModel         = getECfromDatabase(ecModel,noEC);

% Gather kcat values from BRENDA
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% Gather metabolite SMILES, required for DLKcat as kcat source
[ecModel, noSmiles] = findMetSmiles(ecModel);

% Prepare DLKcat input file
writeDLKcatInput(ecModel,[],[],[],[],true);

% Run DLKcat
runDLKcat();

% Load DLKcat output
kcatList_DLKcat = readDLKcatOutput(ecModel);

% Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

% Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);

% Get kcat values across isozymes
ecModel = getKcatAcrossIsozymes(ecModel);

% Get standard kcat
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% Apply kcat constraints from ecModel.ec.kcat to ecModel.S
ecModel = applyKcatConstraints(ecModel);

% Set upper bound of protein pool
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;

saveEcModel(ecModel,'ecKmarx_stage2.yml');

%% STAGE 3: model tuning
% Test maximum growth rate
growthTest = setParam(ecModel,'lb','r_1726',-1000);

% And set growth maximization as the objective function.
growthTest = setParam(growthTest,'obj','r_1912',1);

% Run FBA.
sol = solveLP(growthTest,1);
bioRxnIdx = getIndexes(growthTest,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))

% Relax protein pool constraint
growthTest = setParam(growthTest, 'lb', 'r_1912', 0.48);
growthTest = setParam(growthTest, 'lb', 'prot_pool_exchange', -1000);
growthTest = setParam(growthTest, 'obj', 'prot_pool_exchange', 1);
sol = solveLP(growthTest);

protPoolIdx = strcmp(growthTest.rxns, 'prot_pool_exchange');
fprintf('Protein pool usage is: %.0f mg/gDCW.\n', abs(sol.x(protPoolIdx)))
growthTest = setParam(growthTest,'lb',protPoolIdx,sol.x(protPoolIdx));

% Revert back growth constraint and objective function.
growthTest = setParam(growthTest,'lb','r_1912',0);
growthTest = setParam(growthTest,'obj','r_1912',1);
ecModel = growthTest;

% STEP 43-44 Sensitivity tuning
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

% Inspect the tunedKcats structure in table format.
struct2table(tunedKcats)

saveEcModel(ecModel,'ecKmarx_stage3.yml');

saveEcModel(ecModel,'ecKmarx.yml');