% Generate iDopaNeuro models with different parameters using the
% function XomicsToMultipleModels. 
%
% The batch models are saved in 
% workingDir/multidimensionalModelGeneration


%% Set directory and specificData

clearvars -except workingDir CBTDIR

dataFolder = [CBTDIR filesep 'tutorials' filesep 'dataIntegration' filesep ...
    'XomicsToModel' filesep 'iDopaNeuro' filesep 'data'];

% Omics data filenames
bibliomicData = 'bibliomicData.xlsx';
exometabolomicData = 'exometabolomicData.txt';
transcriptomicData = 'transcriptomicData.txt';

% % % % generic model
% % % inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
% % %     filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
% % % genericModelName = 'Recon3DModel_301_xomics_input.mat';


%% Automatic processing from here on

% Bibliomic data
specificData = preprocessingOmicsModel([dataFolder filesep bibliomicData], 1, 1);

% Exometabolomic data
specificData.exoMet = readtable([dataFolder filesep exometabolomicData]);

% Transcriptomic data
specificData.transcriptomicData = readtable([dataFolder filesep transcriptomicData]);
specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);

% Load generic model
load([CBTDIR filesep 'papers' filesep '2023_iDopaNeuro' filesep 'Recon3DModel_301_thermo.mat'])
modelGenerationConditions.genericModel.model = model;


%action = 'debug';
action = 'batch';
% action = 'iDopaNeuroC';
% action = 'iDopaNeuroCT';

switch action
    case 'batch'
        modelsDir = workingDir;
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end

        % modelGenerationConditions
        modelGenerationConditions.outputDir = modelsDir;        
        % Note: if the data does not vary it is enough to declare them in param.
        % Here the non-varying conditions are left for demonstrative purposes
        modelGenerationConditions.cobraSolver = {'mosek'}; 
        modelGenerationConditions.activeGenesApproach = {'deleteModelGenes', 'oneRxnPerActiveGene'}; 
        modelGenerationConditions.transcriptomicThreshold = [0 2]; % [0 1 2];
        modelGenerationConditions.closeIons = [true false]; 
        modelGenerationConditions.tissueSpecificSolver = {'fastCore', 'thermoKernel'}; 
        modelGenerationConditions.inactiveGenesTranscriptomics = [true false]; 
        modelGenerationConditions.limitBounds = [1e4 1e5]; 
        modelGenerationConditions.curationOverOmics = [true false]; 
        modelGenerationConditions.activeOverInactive = [true false];        
        param.debug = false;
        param.printLevel = 1;
        
    case 'iDopaNeuroC'
        % thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.10000_inactiveGenesT_closedIons_omicsOverCuration
        modelsDir = [workingDir filesep 'iDopaNeuroC'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end
        
        % modelGenerationConditions
        modelGenerationConditions.outputDir = modelsDir;        
        modelGenerationConditions.cobraSolver = {'mosek'};
        modelGenerationConditions.activeGenesApproach = {'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = 2;
        modelGenerationConditions.boundsToRelaxExoMet = {'both'};
        modelGenerationConditions.closeIons = true;
        modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = true; 
        modelGenerationConditions.limitBounds = 10000;
        modelGenerationConditions.curationOverOmics = false; 
        param.debug = true;
        param.printLevel = 2;
        
    case 'iDopaNeuroCT'
        % thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.10000_inactiveGenesT_closedIons_curationOverOmics
        modelsDir = [workingDir filesep 'iDopaNeuroCT'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end

        % modelGenerationConditions
        modelGenerationConditions.outputDir = modelsDir;        
        modelGenerationConditions.cobraSolver = {'ibm_cplex'};
        modelGenerationConditions.activeGenesApproach = {'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = 2;
        modelGenerationConditions.boundsToRelaxExoMet = {'both'};
        modelGenerationConditions.closeIons = true;
        modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = true; 
        modelGenerationConditions.limitBounds = 10000;
        modelGenerationConditions.curationOverOmics = true;
        param.debug = true;
        param.printLevel = 2;

    case 'debug'
        
        modelsDir = [workingDir filesep 'debug'];
        if ~isfolder(modelsDir)
            mkdir(modelsDir);
        end
                            
        %debug generation of a single model type
        % modelGenerationConditions
        % Note: if the data does not vary it is enough to declare them in param.
        % Here the non-varying conditions are left for demonstrative purposes
        modelGenerationConditions.outputDir = modelsDir;
        modelGenerationConditions.cobraSolver = {'gurobi'}; %  {'ibm_cplex', 'gurobi'};
        modelGenerationConditions.activeGenesApproach = {'deleteModelGenes'}; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
        modelGenerationConditions.transcriptomicThreshold = 0; % [0 1 2];
        modelGenerationConditions.boundsToRelaxExoMet = {'both'}; % {'lower', 'both', 'upper'}
        modelGenerationConditions.closeIons = true; % [true false];
        modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'}; % {'fastCore', 'thermoKernel'};
        modelGenerationConditions.inactiveGenesTranscriptomics = false; % [true false];
        modelGenerationConditions.limitBounds = 10000; % [1e3 1e4 1e5]; %was inf
        modelGenerationConditions.curationOverOmics = false; % [true false];
        param.debug = true;
        param.printLevel = 2;
end

% specificData
modelGenerationConditions.specificData.specificData = specificData;

%% Fixed options
param.setObjective = ''; % No objective function
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';
param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true; % Cell culture information
param.findThermoConsistentFluxSubset = 1;

%% Create models

% Remove expressionRxns
if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end

% Create models
directoriesWithModels = XomicsToMultipleModels(modelGenerationConditions, param);

disp('Models saved in:')
disp(modelsDir)
