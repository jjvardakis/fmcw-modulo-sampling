function summary = run_experiments(mode, experiment)
%RUN_EXPERIMENTS Run the report experiments and save their results.
%   run_experiments('quick') runs small execution checks.
%   run_experiments('full') runs all report Monte Carlo trials.
%   run_experiments('full', 'exp2a_iq_imbalance') runs one experiment.
%
% Outputs: results/<mode>/<experiment>/results.mat, console.txt and
% PDF/PNG/FIG files for every figure. A run summary records execution errors.
% An execution pass does not mean that every signal was reconstructed.

if nargin < 1
    mode = 'quick';
end
mode = validatestring(mode, {'quick', 'full'});
experiments = {'exp1_clip_vs_modulo', 'exp2a_iq_imbalance', ...
    'exp2b_noise_floor', 'exp3_snr_sweep', 'exp4_regime_map'};
if nargin >= 2
    experiment = validatestring(experiment, experiments);
    experiments = {experiment};
end

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
outputRoot = fullfile(projectRoot, 'results', mode);
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end
oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restoreVisibility = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVisibility)); %#ok<NASGU>

completed = false(numel(experiments), 1);
seconds = zeros(numel(experiments), 1);
message = strings(numel(experiments), 1);
for index = 1:numel(experiments)
    name = experiments{index};
    folder = fullfile(outputRoot, name);
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    fprintf('Starting %s (%s)\n', name, mode);
    startTime = tic;
    try
        execute_one(fullfile(projectRoot, [name '.m']), folder, strcmp(mode, 'quick'));
        completed(index) = true;
        message(index) = "completed";
    catch err
        message(index) = string(err.message);
        fid = fopen(fullfile(folder, 'error.txt'), 'w');
        fprintf(fid, '%s', getReport(err, 'extended', 'hyperlinks', 'off'));
        fclose(fid);
    end
    seconds(index) = toc(startTime);
    close all;
    fprintf('Finished %s: completed=%d, %.1f seconds\n', name, completed(index), seconds(index));
end
summary = table(string(experiments(:)), completed, seconds, message, ...
    'VariableNames', {'experiment', 'completed', 'seconds', 'message'});
if nargin >= 2
    summaryName = [experiment '_run_summary.csv'];
else
    summaryName = 'run_summary.csv';
end
writetable(summary, fullfile(outputRoot, summaryName));
disp(summary);
if ~all(completed)
    warning('run_experiments:Incomplete', 'Some experiments failed; inspect the saved error logs.');
end
end

function execute_one(scriptPath, outputFolder, v2_quick)
% Run each script in a separate function workspace.
diary(fullfile(outputFolder, 'console.txt'));
stopDiary = onCleanup(@() diary('off')); %#ok<NASGU>
run(scriptPath);

% Save data only: omit graphics, handles, path strings and runner internals.
variables = whos;
results = struct;
for variableIndex = 1:numel(variables)
    name = variables(variableIndex).name;
    value = eval(name);
    if (isnumeric(value) || islogical(value) || isstruct(value)) && ~isempty(value)
        results.(name) = value;
    end
end
save(fullfile(outputFolder, 'results.mat'), '-struct', 'results');

figures = findall(groot, 'Type', 'figure');
for figureIndex = 1:numel(figures)
    name = regexprep(get(figures(figureIndex), 'Name'), '[^a-zA-Z0-9]+', '_');
    fileBase = fullfile(outputFolder, name);
    exportgraphics(figures(figureIndex), [fileBase '.pdf'], 'ContentType', 'vector');
    exportgraphics(figures(figureIndex), [fileBase '.png'], 'Resolution', 180);
    savefig(figures(figureIndex), [fileBase '.fig']);
end
end
