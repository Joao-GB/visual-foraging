function [estParams, expectedX, estX] = runFitSimulation( ...
    marginalizeParams, NumTrials, paramsGen, stairLevel, modelParams)
% runFitSimulation
% Plot-free simulation using PAL_AMRF_setupRF / PAL_AMRF_updateRF.
%
% Inputs:
%   marginalizeParams : kept for compatibility with the previous API
%   NumTrials         : number of trials to run
%   paramsGen         : true generating parameters [alpha beta gamma lambda]
%   stairLevel        : target level used by PAL_CumulativeNormal_Level
%   modelParams       : [beta gamma lambda] used by the AMRF model
%
% Outputs:
%   estParams         : [estAlpha modelBeta modelGamma modelLambda]
%   expectedX         : true x at stairLevel for the generating PF
%   estX              : fitted alpha estimate

    if nargin < 1 || isempty(marginalizeParams), marginalizeParams = 4; end %#ok<NASGU>
    if nargin < 2 || isempty(NumTrials), NumTrials = 150; end
    if nargin < 3 || isempty(paramsGen), paramsGen = [-50 0.15 0.5 0.03]; end
    if nargin < 4 || isempty(stairLevel), stairLevel = 0.75; end
    if nargin < 5 || isempty(modelParams), modelParams = [3 0.5 0.01]; end

    if numel(modelParams) ~= 3
        error('modelParams must be [beta gamma lambda].');
    end

    % True generating parameters
%     alphaTrue  = paramsGen(1);
%     betaTrue   = paramsGen(2);
%     gammaTrue  = paramsGen(3);
    lambdaTrue = paramsGen(4);

    % AMRF model parameters
    betaModel   = modelParams(1);
    gammaModel   = modelParams(2);
    lambdaModel  = modelParams(3);

    % True threshold at requested stair level under the generating PF
    expectedX = PAL_CumulativeNormal( ...
        paramsGen, ...
        min(stairLevel, 1 - lambdaTrue - 1e-4), ...
        'inverse');

    % For simulation, make alpha represent the requested threshold level
    paramsGen(1) = expectedX;

    % Random seed
    if exist('RandStream.m','file')
        s = RandStream.create('mt19937ar','seed','shuffle');
        RandStream.setGlobalStream(s);
    end

    % Reparameterized PF: alpha is the requested stairLevel threshold
    PF = @(params, x, varargin) PAL_CumulativeNormal_Level( ...
        params, ...
        x, stairLevel, varargin{:});

    % Stimulus range used by the adaptive method
    stimRange = single(linspace(-80, -15, 10));

    % Prior over alpha only
    priorAlphaRange = single(linspace(-80, -20, 60));
    prior = PAL_pdfNormal(priorAlphaRange, -50, 20);
    prior = prior ./ sum(prior);

    % AMRF setup
    RF = PAL_AMRF_setupRF( ...
        'priorAlphaRange', priorAlphaRange, ...
        'prior', prior, ...
        'PF', PF, ...
        'beta', betaModel, ...
        'gamma', gammaModel, ...
        'lambda', lambdaModel, ...
        'xMin', min(stimRange), ...
        'xMax', max(stimRange), ...
        'meanmode', 'mean', ...
        'stopCriterion', 'trials', ...
        'stopRule', NumTrials);

    % Trial loop
    while RF.stop ~= 1
        response = rand(1) < PF(paramsGen, RF.xCurrent);
        RF = PAL_AMRF_updateRF(RF, RF.xCurrent, response);
    end

    % Extract estimate from AMRF structure
    estX = getEstimateField(RF);

    % Keep output shape similar, but document it as model settings,
    % not estimated beta/gamma/lambda.
    estParams = [estX, betaModel, gammaModel, lambdaModel];

end

function est = getEstimateField(RF)
% Robustly extract the fitted alpha / threshold estimate from AMRF output.

    candidateFields = {'threshold', 'alpha', 'xCurrent'};
    for i = 1:numel(candidateFields)
        f = candidateFields{i};
        if isfield(RF, f) && ~isempty(RF.(f))
            v = RF.(f);
            est = v(end);
            return;
        end
    end

    error('Could not find an estimate field in the RF structure.');
end