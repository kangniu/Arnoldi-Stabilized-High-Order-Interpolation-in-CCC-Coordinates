function results = paper_numerical_experiments_all_in_one()
%PAPER_NUMERICAL_EXPERIMENTS_ALL_IN_ONE Reproduce all numerical examples.
%
% This single MATLAB file contains the complete code used for the numerical
% experiments in the manuscript.  It is a self-contained version of the
% modular code in matlab/experiments and matlab/common.
%
% Run from MATLAB with:
%   results = paper_numerical_experiments_all_in_one;
%
% Outputs are written to the output/ directory next to this file.

close all; clc;

matlabRoot = fileparts(mfilename('fullpath'));
if isempty(matlabRoot)
    matlabRoot = pwd;
end

cfg = default_experiment_config(matlabRoot);
cfg.writePaperFigures = false;  % Set true to overwrite root PNGs used by main.tex.
cfg.makePlots = true;
cfg.saveData = true;

rng(cfg.seed, 'twister');

fprintf('Running all Grassmann interpolation experiments...\n');
fprintf('Output directory: %s\n\n', cfg.outputDir);

resultsExample1 = run_example1_dual_stability(cfg);
resultsExample2 = run_example2_noise(cfg);
resultsExample3 = run_example3_helmholtz_rom(cfg);

if cfg.saveData
    save(fullfile(cfg.outputDataDir, 'all_experiment_results.mat'), ...
        'cfg', 'resultsExample1', 'resultsExample2', 'resultsExample3');
end

results = struct();
results.cfg = cfg;
results.example1 = resultsExample1;
results.example2 = resultsExample2;
results.example3 = resultsExample3;

fprintf('\nFinished all experiments.\n');

if nargout == 0
    clear results;
end
end

%% Source: .\matlab\experiments\run_example1_dual_stability.m

function out = run_example1_dual_stability(cfg)
%RUN_EXAMPLE1_DUAL_STABILITY High-order stability tests on Gr(1000,10).

fprintf('\nExample 1: high-order interpolation and dual stability\n');

n = 1000;
p = 10;
curve = synthetic_curve_setup(n, p, cfg.syntheticSeed);
tEval = linspace(0, 1, cfg.evalPoints).';
trueSample = synthetic_curve_sample(curve, tEval);

methods = make_method_list('full');
paperMethodNames = {'normal_direct', 'raw_direct', 'maxvol_direct', 'ccc_cva'};

out = struct();
out.degree15 = run_one_case(cfg, curve, trueSample, tEval, 8, 4, ...
    'example1_degree15', 'dual_stability_figure.png', methods, paperMethodNames);
out.degree35 = run_one_case(cfg, curve, trueSample, tEval, 18, [], ...
    'example1_degree35', 'dual_stability_figure2.png', methods, paperMethodNames);
end

function caseOut = run_one_case(cfg, curve, trueSample, tEval, m, centerIndex, prefix, paperFigureName, methods, paperMethodNames)
tTrain = cheb_lobatto_nodes(m);
if isempty(centerIndex)
    centerIndex = closest_center_index(tTrain, 0.5);
end
trainSample = synthetic_curve_sample(curve, tTrain);

fprintf('  Nodes m = %d, Hermite degree = %d, center t = %.4f\n', ...
    m, 2*m - 1, tTrain(centerIndex));

comparison = run_hermite_comparison(tTrain, trainSample, tEval, trueSample, centerIndex, methods);
save_comparison_outputs(cfg, prefix, comparison);

chartTable = chart_condition_table(tTrain, trainSample.U, centerIndex);
writetable(chartTable, fullfile(cfg.outputTableDir, [prefix, '_chart_condition.csv']));

if cfg.makePlots
    fullFig = fullfile(cfg.outputFigureDir, [prefix, '_full.png']);
    plot_error_orthogonality(tEval, comparison.results, ...
        sprintf('Example 1, m=%d, Hermite degree %d: full comparison', m, 2*m-1), fullFig);

    paperFig = fullfile(cfg.outputFigureDir, paperFigureName);
    paperResults = select_results(comparison.results, paperMethodNames);
    plot_error_orthogonality(tEval, paperResults, ...
        sprintf('Example 1, m=%d, Hermite degree %d', m, 2*m-1), paperFig);
    write_paper_figure_if_requested(cfg, paperFig, paperFigureName);
end

caseOut = struct();
caseOut.tTrain = tTrain;
caseOut.centerIndex = centerIndex;
caseOut.comparison = comparison;
caseOut.chartCondition = chartTable;
end


%% Source: .\matlab\experiments\run_example2_noise.m

function out = run_example2_noise(cfg)
%RUN_EXAMPLE2_NOISE Sensitivity to perturbed input data.

fprintf('\nExample 2: sensitivity to data perturbations\n');

n = 1000;
p = 10;
m = 10;
noiseLevel = 1e-10;

curve = synthetic_curve_setup(n, p, cfg.syntheticSeed);
tTrain = cheb_lobatto_nodes(m);
tEval = linspace(0, 1, cfg.evalPoints).';
centerIndex = closest_center_index(tTrain, 0.5);

trainSample = synthetic_curve_sample(curve, tTrain, ...
    struct('noiseLevel', noiseLevel, 'noiseSeed', cfg.noiseSeed));
trueSample = synthetic_curve_sample(curve, tEval);

fprintf('  Nodes m = %d, Hermite degree = %d, center t = %.4f, noise = %.1e\n', ...
    m, 2*m - 1, tTrain(centerIndex), noiseLevel);

methods = make_method_list('full');
paperMethodNames = {'normal_direct', 'raw_direct', 'maxvol_direct', 'ccc_cva'};
comparison = run_hermite_comparison(tTrain, trainSample, tEval, trueSample, centerIndex, methods);
comparison.noiseLevel = noiseLevel;

prefix = 'example2_noise';
save_comparison_outputs(cfg, prefix, comparison);

if cfg.makePlots
    fullFig = fullfile(cfg.outputFigureDir, [prefix, '_full.png']);
    plot_error_only(tEval, comparison.results, ...
        sprintf('Example 2, noise %.1e: full comparison', noiseLevel), fullFig, 't');

    paperFig = fullfile(cfg.outputFigureDir, 'noise_sensitivity_figure.png');
    paperResults = select_results(comparison.results, paperMethodNames);
    plot_error_only(tEval, paperResults, ...
        sprintf('Example 2, noise %.1e', noiseLevel), paperFig, 't');
    write_paper_figure_if_requested(cfg, paperFig, 'noise_sensitivity_figure.png');
end

out = struct();
out.tTrain = tTrain;
out.centerIndex = centerIndex;
out.noiseLevel = noiseLevel;
out.comparison = comparison;
end


%% Source: .\matlab\experiments\run_example3_helmholtz_rom.m

function out = run_example3_helmholtz_rom(cfg)
%RUN_EXAMPLE3_HELMHOLTZ_ROM Parametric 1D Helmholtz ROM interpolation.

fprintf('\nExample 3: parametric 1D Helmholtz ROM\n');

n = 500;
p = 8;
m = 12;

curve = helmholtz_curve_setup(n, p, cfg.helmholtzSeed);
tTrain = cheb_lobatto_nodes(m);
tEval = linspace(0, 1, cfg.evalPoints).';
centerIndex = closest_center_index(tTrain, 0.5);

trainSample = helmholtz_curve_sample(curve, tTrain);
trueSample = helmholtz_curve_sample(curve, tEval);

fprintf('  Nodes m = %d, Hermite degree = %d, center k = %.4f\n', ...
    m, 2*m - 1, trainSample.k(centerIndex));

methods = make_method_list('full');
paperMethodNames = {'normal_direct', 'raw_direct', 'maxvol_direct', 'ccc_cva'};
comparison = run_hermite_comparison(tTrain, trainSample, tEval, trueSample, centerIndex, methods);
comparison.kEval = trueSample.k;
comparison.kTrain = trainSample.k;

prefix = 'example3_helmholtz_rom';
save_comparison_outputs(cfg, prefix, comparison);

if cfg.makePlots
    fullFig = fullfile(cfg.outputFigureDir, [prefix, '_full.png']);
    plot_error_only(trueSample.k, comparison.results, ...
        'Example 3, Helmholtz ROM: full comparison', fullFig, 'k');

    paperFig = fullfile(cfg.outputFigureDir, 'helmholtz_rom_error.png');
    paperResults = select_results(comparison.results, paperMethodNames);
    plot_error_only(trueSample.k, paperResults, ...
        'Example 3, Helmholtz ROM', paperFig, 'k');
    write_paper_figure_if_requested(cfg, paperFig, 'helmholtz_rom_error.png');
end

out = struct();
out.tTrain = tTrain;
out.kTrain = trainSample.k;
out.centerIndex = centerIndex;
out.comparison = comparison;
end


%% Source: .\matlab\common\apply_chart_inverse.m

function U = apply_chart_inverse(chart, V)
%APPLY_CHART_INVERSE Map a chart-adapted representative back to original rows.

switch chart.type
    case {'raw', 'maxvol'}
        U = zeros(chart.n, chart.p);
        U(chart.perm, :) = V;
    case 'ccc'
        U = chart.Q * V;
    otherwise
        error('Unknown chart type "%s".', chart.type);
end
end


%% Source: .\matlab\common\apply_chart_transform.m

function V = apply_chart_transform(chart, U)
%APPLY_CHART_TRANSFORM Transform a Stiefel/tangent matrix into chart rows.

switch chart.type
    case {'raw', 'maxvol'}
        V = U(chart.perm, :);
    case 'ccc'
        V = chart.Q.' * U;
    otherwise
        error('Unknown chart type "%s".', chart.type);
end
end


%% Source: .\matlab\common\build_chart.m

function chart = build_chart(chartType, Uref)
%BUILD_CHART Construct a fixed Grassmann affine chart.
%
% chartType can be:
%   'raw'    first p rows are used as the affine block;
%   'maxvol' greedy maximum-volume row selection;
%   'ccc'    canonical centering by a complete QR transform of Uref.

chartType = lower(string(chartType));
[n, p] = size(Uref);

chart = struct();
chart.type = char(chartType);
chart.n = n;
chart.p = p;

switch chartType
    case "raw"
        chart.perm = 1:n;
        chart.Q = [];
    case "maxvol"
        selected = maxvol_rows(Uref);
        rest = setdiff(1:n, selected, 'stable');
        chart.perm = [selected(:).', rest(:).'];
        chart.selectedRows = selected(:);
        chart.Q = [];
    case "ccc"
        [Q, ~] = qr(Uref);
        chart.Q = Q;
        chart.perm = [];
    otherwise
        error('Unknown chart type "%s".', chartType);
end
end


%% Source: .\matlab\common\chart_condition_table.m

function T = chart_condition_table(t, Ucell, centerIndex)
%CHART_CONDITION_TABLE Frobenius chart condition for raw, maxvol, and CCC.

Uref = Ucell{centerIndex};
charts = {build_chart('raw', Uref), build_chart('maxvol', Uref), build_chart('ccc', Uref)};
names = {'RawFixedBlock', 'MaxvolRowSelection', 'CCCCentered'};

data = zeros(numel(t), numel(charts));
for j = 1:numel(charts)
    coord = map_chart_data(charts{j}, Ucell, zero_tangent_cells(Ucell));
    data(:, j) = coord.kappaF;
end

T = table(t(:), data(:, 1), data(:, 2), data(:, 3), ...
    'VariableNames', [{'t'}, names]);
end


%% Source: .\matlab\common\cheb_lobatto_nodes.m

function t = cheb_lobatto_nodes(m)
%CHEB_LOBATTO_NODES Chebyshev-Lobatto nodes mapped to [0,1].

if m < 2
    t = 0.5;
    return;
end

j = (0:m-1).';
t = 0.5 - 0.5*cos(pi*j/(m-1));
end


%% Source: .\matlab\common\closest_center_index.m

function idx = closest_center_index(t, target)
%CLOSEST_CENTER_INDEX Index of the node closest to target.

if nargin < 2
    target = 0.5;
end
[~, idx] = min(abs(t(:) - target));
end


%% Source: .\matlab\common\confluent_vandermonde.m

function C = confluent_vandermonde(t, degree)
%CONFLUENT_VANDERMONDE Hermite matrix for values followed by derivatives.

t = t(:);
m = numel(t);
C = zeros(2*m, degree + 1);
C(1:m, :) = vandermonde_values(t, degree);

for j = 2:degree+1
    C(m+1:end, j) = (j - 1) * t.^(j - 2);
end
end


%% Source: .\matlab\common\cva_arnoldi.m

function arn = cva_arnoldi(t, degree)
%CVA_ARNOLDI Scaled Arnoldi basis for confluent Hermite interpolation.

t = t(:);
m = numel(t);
N = 2*m;
if degree > N - 1
    error('CVA degree %d exceeds the Hermite data dimension %d.', degree, N - 1);
end

Q = zeros(N, degree + 1);
H = zeros(degree + 1, degree);
Q(:, 1) = [ones(m, 1); zeros(m, 1)];

for j = 1:degree
    qTop = Q(1:m, j);
    qBot = Q(m+1:end, j);
    v = [t.*qTop; qTop + t.*qBot];

    for pass = 1:2
        for i = 1:j
            hij = (Q(:, i).'*v) / m;
            H(i, j) = H(i, j) + hij;
            v = v - Q(:, i)*hij;
        end
    end

    hnext = norm(v) / sqrt(m);
    H(j + 1, j) = hnext;
    if hnext <= 100*eps
        warning('cva_arnoldi:Breakdown', ...
            'Near breakdown at step %d; remaining basis columns are set to NaN.', j);
        Q(:, j+1:end) = NaN;
        break;
    end
    Q(:, j + 1) = v / hnext;
end

arn = struct();
arn.Q = Q;
arn.H = H;
arn.degree = degree;
arn.m = m;
arn.etaQ = norm((Q.'*Q)/m - eye(degree + 1), 2);
arn.condQ = cond(Q);
arn.minSubdiag = min(abs(diag(H, -1)));
end


%% Source: .\matlab\common\default_experiment_config.m

function cfg = default_experiment_config(matlabRoot)
%DEFAULT_EXPERIMENT_CONFIG Shared settings and output locations.

if nargin < 1 || isempty(matlabRoot)
    matlabRoot = fileparts(fileparts(mfilename('fullpath')));
end

cfg = struct();
cfg.seed = 20260716;
cfg.syntheticSeed = 1137;
cfg.noiseSeed = 2718;
cfg.helmholtzSeed = 3141;
cfg.evalPoints = 200;
cfg.makePlots = true;
cfg.saveData = true;
cfg.writePaperFigures = false;

cfg.matlabRoot = matlabRoot;
cfg.projectRoot = fileparts(matlabRoot);
cfg.outputDir = fullfile(matlabRoot, 'output');
cfg.outputFigureDir = fullfile(cfg.outputDir, 'figures');
cfg.outputTableDir = fullfile(cfg.outputDir, 'tables');
cfg.outputDataDir = fullfile(cfg.outputDir, 'data');

ensure_dir(cfg.outputDir);
ensure_dir(cfg.outputFigureDir);
ensure_dir(cfg.outputTableDir);
ensure_dir(cfg.outputDataDir);
end


%% Source: .\matlab\common\empty_alg.m

function alg = empty_alg(typeName, degree)
%EMPTY_ALG Placeholder algorithm diagnostics after a failed method.

alg = struct();
alg.type = typeName;
alg.degree = degree;
alg.condSystem = NaN;
alg.coeffNormF = NaN;
alg.etaQ = NaN;
alg.condQ = NaN;
alg.minSubdiag = NaN;
alg.basisGrowth = NaN;
alg.recurrenceDiagnostic = NaN;
end


%% Source: .\matlab\common\ensure_dir.m

function ensure_dir(pathName)
%ENSURE_DIR Create a directory if it does not already exist.

if ~exist(pathName, 'dir')
    mkdir(pathName);
end
end


%% Source: .\matlab\common\eval_arnoldi_basis.m

function Phi = eval_arnoldi_basis(s, H)
%EVAL_ARNOLDI_BASIS Evaluate a V+A scalar basis from Hessenberg coefficients.

s = s(:);
M = numel(s);
degree = size(H, 2);
Phi = zeros(M, degree + 1);
Phi(:, 1) = ones(M, 1);

for j = 1:degree
    v = s.*Phi(:, j) - Phi(:, 1:j)*H(1:j, j);
    hnext = H(j + 1, j);
    if abs(hnext) <= realmin
        Phi(:, j+1:end) = NaN;
        break;
    end
    Phi(:, j + 1) = v / hnext;
end
end


%% Source: .\matlab\common\eval_augmented_basis.m

function [PhiValue, PhiDeriv, PhiAug] = eval_augmented_basis(s, H)
%EVAL_AUGMENTED_BASIS Evaluate CV+A basis and derivative basis at points s.

s = s(:);
M = numel(s);
degree = size(H, 2);
PhiAug = zeros(2*M, degree + 1);
PhiAug(:, 1) = [ones(M, 1); zeros(M, 1)];

for j = 1:degree
    top = PhiAug(1:M, j);
    bot = PhiAug(M+1:end, j);
    v = [s.*top; top + s.*bot] - PhiAug(:, 1:j)*H(1:j, j);
    hnext = H(j + 1, j);
    if abs(hnext) <= realmin
        PhiAug(:, j+1:end) = NaN;
        break;
    end
    PhiAug(:, j + 1) = v / hnext;
end

PhiValue = PhiAug(1:M, :);
PhiDeriv = PhiAug(M+1:end, :);
end


%% Source: .\matlab\common\evaluate_hermite_method.m

function result = evaluate_hermite_method(method, tTrain, Utrain, Dtrain, tEval, Utrue, centerIndex)
%EVALUATE_HERMITE_METHOD Fit one Hermite method and evaluate projector error.

if nargin < 7 || isempty(centerIndex)
    centerIndex = ceil(numel(tTrain)/2);
end

tTrain = tTrain(:);
tEval = tEval(:);
Uref = Utrain{centerIndex};
[n, p] = size(Uref);

result = struct();
result.name = method.name;
result.label = method.label;
result.coordType = method.coordType;
result.algebra = method.algebra;
result.tEval = tEval;
result.error = NaN(numel(tEval), 1);
result.orthLoss = NaN(numel(tEval), 1);
result.failure = '';

try
    switch lower(method.coordType)
        case 'normal'
            coord = map_normal_data(Uref, Utrain, Dtrain);
            reconstruct = @(xrow) reconstruct_normal_point(coord, xrow);
        case {'raw', 'maxvol', 'ccc'}
            chart = build_chart(method.coordType, Uref);
            coord = map_chart_data(chart, Utrain, Dtrain);
            reconstruct = @(xrow) reconstruct_chart_point(chart, xrow);
        otherwise
            error('Unknown coordinate type "%s".', method.coordType);
    end

    switch lower(method.algebra)
        case 'direct'
            [Xeval, alg] = hermite_interpolate_direct(tTrain, coord.X, coord.Xdot, tEval);
        case 'cva'
            [Xeval, alg] = hermite_interpolate_cva(tTrain, coord.X, coord.Xdot, tEval);
        otherwise
            error('Unknown algebra type "%s".', method.algebra);
    end

    for j = 1:numel(tEval)
        Uhat = reconstruct(Xeval(j, :));
        result.error(j) = projector_relative_error(Uhat, Utrue{j});
        result.orthLoss(j) = orthogonality_loss(Uhat);
    end

    result.coord = strip_large_coord(coord);
    result.alg = alg;
catch ME
    result.failure = ME.message;
    result.coord = struct();
    result.alg = empty_alg(method.algebra, 2*numel(tTrain)-1);
    warning('evaluate_hermite_method:Failed', ...
        'Method %s failed: %s', method.name, ME.message);
end

result.summary = summarize_method_result(result, n, p, numel(tTrain));
end


%% Source: .\matlab\common\floor_log_values.m

function y = floor_log_values(x)
%FLOOR_LOG_VALUES Avoid zeros and negative roundoff in logarithmic plots.

y = x(:);
y(~isfinite(y)) = NaN;
y = max(y, 1e-16);
end


%% Source: .\matlab\common\grassmann_exp.m

function U = grassmann_exp(U0, X)
%GRASSMANN_EXP Grassmann exponential from a horizontal tangent matrix.

X = X - U0*(U0.'*X);
[A, S, B] = svd(X, 'econ');
sigma = diag(S);
C = diag(cos(sigma));
Sn = diag(sin(sigma));
U = U0*(B*C*B.') + A*Sn*B.';
end


%% Source: .\matlab\common\grassmann_log.m

function X = grassmann_log(U0, U)
%GRASSMANN_LOG Grassmann logarithm in the horizontal tangent at U0.

M = (U - U0*(U0.'*U)) / (U0.'*U);
[A, S, B] = svd(M, 'econ');
theta = atan(diag(S));
X = A*diag(theta)*B.';
X = X - U0*(U0.'*X);
end


%% Source: .\matlab\common\helmholtz_curve_sample.m

function sample = helmholtz_curve_sample(curve, t)
%HELMHOLTZ_CURVE_SAMPLE Stiefel data for k in [10,20] normalized by t in [0,1].

t = t(:);
m = numel(t);
U = cell(m, 1);
Delta = cell(m, 1);
Y = cell(m, 1);
Ydot = cell(m, 1);
k = curve.kMin + (curve.kMax - curve.kMin)*t;
dkdt = curve.kMax - curve.kMin;
I = speye(curve.n);

for i = 1:m
    A = curve.L - (k(i)^2)*I;
    Yi = A \ curve.F;
    YiPrimeK = A \ (2*k(i)*Yi);
    YiPrimeT = dkdt*YiPrimeK;
    [Ui, Di] = horizontal_lift_from_full_rank(Yi, YiPrimeT);

    U{i} = Ui;
    Delta{i} = Di;
    Y{i} = Yi;
    Ydot{i} = YiPrimeT;
end

sample = struct();
sample.t = t;
sample.k = k;
sample.U = U;
sample.Delta = Delta;
sample.Y = Y;
sample.Ydot = Ydot;
sample.n = curve.n;
sample.p = curve.p;
end


%% Source: .\matlab\common\helmholtz_curve_setup.m

function curve = helmholtz_curve_setup(n, p, seed)
%HELMHOLTZ_CURVE_SETUP One-dimensional Helmholtz ROM data generator.

if nargin < 3
    seed = 3141;
end

h = 1/(n + 1);
e = ones(n, 1);
L = spdiags([-e, 2*e, -e], [-1, 0, 1], n, n) / h^2;
x = (1:n).' * h;

centers = linspace(0.15, 0.85, p);
width = 0.035;
F = zeros(n, p);
for j = 1:p
    F(:, j) = exp(-((x - centers(j))/width).^2);
    F(:, j) = F(:, j) / max(norm(F(:, j)), realmin);
end

state = rng;
rng(seed, 'twister');
F = F + 0.02*randn(n, p);
for j = 1:p
    F(:, j) = F(:, j) / max(norm(F(:, j)), realmin);
end
rng(state);

curve = struct();
curve.n = n;
curve.p = p;
curve.seed = seed;
curve.h = h;
curve.x = x;
curve.L = L;
curve.F = F;
curve.kMin = 10;
curve.kMax = 20;
end


%% Source: .\matlab\common\hermite_interpolate_cva.m

function [Xeval, alg] = hermite_interpolate_cva(t, X, Xdot, s)
%HERMITE_INTERPOLATE_CVA Confluent V+A interpolation for many RHS columns.

t = t(:);
s = s(:);
degree = 2*numel(t) - 1;
arn = cva_arnoldi(t, degree);
Y = [X; Xdot];

A = arn.Q \ Y;
[PhiValue, ~, PhiAug] = eval_augmented_basis(s, arn.H);
Xeval = PhiValue * A;

alg = struct();
alg.type = 'cva';
alg.degree = degree;
alg.condSystem = NaN;
alg.coeffNormF = norm(A, 'fro');
alg.etaQ = arn.etaQ;
alg.condQ = arn.condQ;
alg.minSubdiag = arn.minSubdiag;
alg.basisGrowth = max(vecnorm(PhiValue, 2, 2));
alg.recurrenceDiagnostic = recurrence_diagnostic(PhiAug, arn.H);
alg.H = arn.H;
end


%% Source: .\matlab\common\hermite_interpolate_direct.m

function [Xeval, alg] = hermite_interpolate_direct(t, X, Xdot, s)
%HERMITE_INTERPOLATE_DIRECT Direct monomial confluent Vandermonde solve.

t = t(:);
s = s(:);
degree = 2*numel(t) - 1;
C = confluent_vandermonde(t, degree);
Y = [X; Xdot];

coef = C \ Y;
V = vandermonde_values(s, degree);
Xeval = V * coef;

alg = struct();
alg.type = 'direct';
alg.degree = degree;
alg.condSystem = cond(C);
alg.coeffNormF = norm(coef, 'fro');
alg.etaQ = NaN;
alg.condQ = NaN;
alg.minSubdiag = NaN;
alg.basisGrowth = max(vecnorm(V, 2, 2));
alg.recurrenceDiagnostic = NaN;
end


%% Source: .\matlab\common\hide_axes_toolbars.m

function hide_axes_toolbars(fig)
%HIDE_AXES_TOOLBARS Prevent MATLAB axes toolbars from appearing in exports.

axesHandles = findall(fig, 'Type', 'axes');
for k = 1:numel(axesHandles)
    ax = axesHandles(k);
    if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end
end
end


%% Source: .\matlab\common\horizontal_lift_from_full_rank.m

function [U, Delta, R] = horizontal_lift_from_full_rank(Y, Ydot)
%HORIZONTAL_LIFT_FROM_FULL_RANK Convert full-rank matrix data to Stiefel data.
%
% If Y = U R is a thin QR factorization, the Grassmann tangent represented
% by Ydot has horizontal lift Delta = (I-UU') Ydot R^{-1}.

[U, R] = thin_qr_positive(Y);
Delta = (Ydot - U*(U.'*Ydot)) / R;
Delta = Delta - U*(U.'*Delta);
end


%% Source: .\matlab\common\lagrange_interpolate_direct.m

function [Xeval, alg] = lagrange_interpolate_direct(t, X, s)
%LAGRANGE_INTERPOLATE_DIRECT Direct monomial Vandermonde interpolation.

t = t(:);
s = s(:);
degree = numel(t) - 1;
Vtrain = vandermonde_values(t, degree);
coef = Vtrain \ X;
Veval = vandermonde_values(s, degree);
Xeval = Veval * coef;

alg = struct();
alg.type = 'direct-lagrange';
alg.degree = degree;
alg.condSystem = cond(Vtrain);
alg.coeffNormF = norm(coef, 'fro');
alg.etaQ = NaN;
alg.condQ = NaN;
alg.minSubdiag = NaN;
alg.basisGrowth = max(vecnorm(Veval, 2, 2));
alg.recurrenceDiagnostic = NaN;
end


%% Source: .\matlab\common\lagrange_interpolate_va.m

function [Xeval, alg] = lagrange_interpolate_va(t, X, s)
%LAGRANGE_INTERPOLATE_VA Vandermonde-with-Arnoldi interpolation.

t = t(:);
s = s(:);
degree = numel(t) - 1;
arn = va_arnoldi(t, degree);
A = arn.Q \ X;
Phi = eval_arnoldi_basis(s, arn.H);
Xeval = Phi * A;

alg = struct();
alg.type = 'va';
alg.degree = degree;
alg.condSystem = NaN;
alg.coeffNormF = norm(A, 'fro');
alg.etaQ = arn.etaQ;
alg.condQ = arn.condQ;
alg.minSubdiag = arn.minSubdiag;
alg.basisGrowth = max(vecnorm(Phi, 2, 2));
alg.recurrenceDiagnostic = NaN;
alg.H = arn.H;
end


%% Source: .\matlab\common\make_method_list.m

function methods = make_method_list(kind)
%MAKE_METHOD_LIST Methods used in the numerical comparisons.

if nargin < 1
    kind = 'full';
end
kind = lower(char(kind));

allMethods = {
    'normal_direct', 'Normal + direct Hermite', 'normal', 'direct', true;
    'raw_direct',    'Raw fixed block + direct Hermite', 'raw', 'direct', true;
    'maxvol_direct', 'Maxvol row selection + direct Hermite', 'maxvol', 'direct', true;
    'ccc_direct',    'CCC + direct Hermite', 'ccc', 'direct', false;
    'normal_cva',    'Normal + CV+A', 'normal', 'cva', false;
    'raw_cva',       'Raw fixed block + CV+A', 'raw', 'cva', false;
    'maxvol_cva',    'Maxvol row selection + CV+A', 'maxvol', 'cva', false;
    'ccc_cva',       'CCC-CV+A', 'ccc', 'cva', true;
    };

switch kind
    case 'paper'
        keep = cell2mat(allMethods(:, 5));
    case 'full'
        keep = true(size(allMethods, 1), 1);
    otherwise
        error('Unknown method-list kind "%s".', kind);
end

rows = allMethods(keep, :);
methods = struct('name', {}, 'label', {}, 'coordType', {}, 'algebra', {});
for i = 1:size(rows, 1)
    methods(i).name = rows{i, 1};
    methods(i).label = rows{i, 2};
    methods(i).coordType = rows{i, 3};
    methods(i).algebra = rows{i, 4};
end
end


%% Source: .\matlab\common\map_chart_data.m

function coord = map_chart_data(chart, Ucell, Dcell)
%MAP_CHART_DATA Map Stiefel samples and horizontal lifts to affine coordinates.

m = numel(Ucell);
n = chart.n;
p = chart.p;
q = n - p;

X = zeros(m, q*p);
Xdot = zeros(m, q*p);
kappaF = zeros(m, 1);
kappa2 = zeros(m, 1);
sigmaMin = zeros(m, 1);
coordNorm = zeros(m, 1);

for i = 1:m
    V = apply_chart_transform(chart, Ucell{i});
    W = apply_chart_transform(chart, Dcell{i});

    V1 = V(1:p, :);
    V2 = V(p+1:end, :);
    W1 = W(1:p, :);
    W2 = W(p+1:end, :);

    Xi = V2 / V1;
    Xidot = (W2 - Xi*W1) / V1;

    X(i, :) = Xi(:).';
    Xdot(i, :) = Xidot(:).';

    invV1 = V1 \ eye(p);
    s = svd(V1);
    sigmaMin(i) = min(s);
    kappaF(i) = norm(invV1, 'fro');
    kappa2(i) = norm(invV1, 2);
    coordNorm(i) = norm(Xi, 'fro');
end

coord = struct();
coord.kind = 'chart';
coord.chart = chart;
coord.X = X;
coord.Xdot = Xdot;
coord.kappaF = kappaF;
coord.kappa2 = kappa2;
coord.sigmaMin = sigmaMin;
coord.coordNorm = coordNorm;
end


%% Source: .\matlab\common\map_normal_data.m

function coord = map_normal_data(Uref, Ucell, Dcell)
%MAP_NORMAL_DATA Map samples into normal coordinates at Uref.

m = numel(Ucell);
[n, p] = size(Uref);
X = zeros(m, n*p);
Xdot = zeros(m, n*p);
coordNorm = zeros(m, 1);
h = 1e-6;

for i = 1:m
    Xi = grassmann_log(Uref, Ucell{i});

    Up = orthonormalize(Ucell{i} + h*Dcell{i});
    Um = orthonormalize(Ucell{i} - h*Dcell{i});
    Xip = grassmann_log(Uref, Up);
    Xim = grassmann_log(Uref, Um);
    Xidot = (Xip - Xim) / (2*h);

    X(i, :) = Xi(:).';
    Xdot(i, :) = Xidot(:).';
    coordNorm(i) = norm(Xi, 'fro');
end

coord = struct();
coord.kind = 'normal';
coord.Uref = Uref;
coord.n = n;
coord.p = p;
coord.X = X;
coord.Xdot = Xdot;
coord.coordNorm = coordNorm;
end


%% Source: .\matlab\common\maxvol_rows.m

function selected = maxvol_rows(A, tol, maxIter)
%MAXVOL_ROWS Greedy maximum-volume row subset for a tall full-rank matrix.
%
% The routine starts from QR-pivoted rows and performs the standard maxvol
% row-exchange iteration until all interpolation coefficients have magnitude
% at most 1+tol.  It is intended for the small p used in these experiments.

if nargin < 2 || isempty(tol)
    tol = 1e-2;
end
if nargin < 3 || isempty(maxIter)
    maxIter = 500;
end

[n, p] = size(A);
if n < p
    error('maxvol_rows:Shape', 'The input matrix must have at least p rows.');
end

[~, ~, piv] = qr(A.', 'vector');
selected = piv(1:p);

for iter = 1:maxIter
    B = A(selected, :);
    if rcond(B) < eps
        break;
    end
    C = A / B;
    [maxVal, linearIndex] = max(abs(C(:)));
    if maxVal <= 1 + tol
        return;
    end
    [rowIndex, colIndex] = ind2sub(size(C), linearIndex);
    selected(colIndex) = rowIndex;
    selected = repair_duplicate_rows(A, selected);
end

warning('maxvol_rows:NoConvergence', ...
    'Maxvol row selection stopped after %d iterations.', maxIter);
end


%% Source: .\matlab\common\method_plot_style.m

function style = method_plot_style(methodName)
%METHOD_PLOT_STYLE Fixed high-contrast plotting style by method name.

style = struct();
style.Color = [0.15, 0.15, 0.15];
style.LineStyle = '-';
style.LineWidth = 1.7;

switch lower(char(methodName))
    case 'normal_direct'
        style.Color = [0.000, 0.250, 0.750];
        style.LineStyle = '--';
        style.LineWidth = 1.8;
    case 'raw_direct'
        style.Color = [0.800, 0.120, 0.120];
        style.LineStyle = '-';
        style.LineWidth = 1.9;
    case 'maxvol_direct'
        style.Color = [0.000, 0.500, 0.180];
        style.LineStyle = '-.';
        style.LineWidth = 1.9;
    case 'ccc_cva'
        style.Color = [0.000, 0.000, 0.000];
        style.LineStyle = '-';
        style.LineWidth = 2.2;
    case 'ccc_direct'
        style.Color = [0.450, 0.450, 0.450];
        style.LineStyle = ':';
        style.LineWidth = 1.8;
    case 'normal_cva'
        style.Color = [0.000, 0.620, 0.720];
        style.LineStyle = '--';
        style.LineWidth = 1.7;
    case 'raw_cva'
        style.Color = [0.950, 0.550, 0.150];
        style.LineStyle = ':';
        style.LineWidth = 1.8;
    case 'maxvol_cva'
        style.Color = [0.350, 0.700, 0.250];
        style.LineStyle = ':';
        style.LineWidth = 1.8;
end
end


%% Source: .\matlab\common\orthogonality_loss.m

function loss = orthogonality_loss(U)
%ORTHOGONALITY_LOSS Frobenius Stiefel orthogonality residual.

if any(~isfinite(U(:)))
    loss = NaN;
    return;
end

p = size(U, 2);
loss = norm(U.'*U - eye(p), 'fro');
end


%% Source: .\matlab\common\orthonormalize.m

function U = orthonormalize(Y)
%ORTHONORMALIZE Return a deterministic thin orthonormal basis.

[U, ~] = thin_qr_positive(Y);
end


%% Source: .\matlab\common\plot_error_only.m

function plot_error_only(tEval, results, titleText, outFile, xLabelText)
%PLOT_ERROR_ONLY Single-panel semilog plot of relative projector error.

if nargin < 5 || isempty(xLabelText)
    xLabelText = 't';
end

fig = figure('Color', 'w', 'Position', [100, 100, 760, 460]);
hold on;
for i = 1:numel(results)
    style = method_plot_style(results(i).name);
    semilogy(tEval, floor_log_values(results(i).error), ...
        'Color', style.Color, 'LineStyle', style.LineStyle, ...
        'LineWidth', style.LineWidth);
end
grid on;
set(gca, 'YScale', 'log');
xlabel(xLabelText);
ylabel('relative projector error');
title(titleText, 'Interpreter', 'none');
legend({results.label}, 'Interpreter', 'none', 'Location', 'southoutside', 'NumColumns', 2);

if nargin >= 4 && ~isempty(outFile)
    hide_axes_toolbars(fig);
    exportgraphics(fig, outFile, 'Resolution', 300);
end
end


%% Source: .\matlab\common\plot_error_orthogonality.m

function plot_error_orthogonality(tEval, results, titleText, outFile)
%PLOT_ERROR_ORTHOGONALITY Two-panel plot for reconstruction and orthogonality.

fig = figure('Color', 'w', 'Position', [100, 100, 1120, 430]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, titleText, 'Interpreter', 'none');

nexttile;
hold on;
for i = 1:numel(results)
    style = method_plot_style(results(i).name);
    semilogy(tEval, floor_log_values(results(i).error), ...
        'Color', style.Color, 'LineStyle', style.LineStyle, ...
        'LineWidth', style.LineWidth);
end
grid on;
set(gca, 'YScale', 'log');
xlabel('t');
ylabel('relative projector error');
legend({results.label}, 'Interpreter', 'none', 'Location', 'southoutside', 'NumColumns', 2);

nexttile;
hold on;
for i = 1:numel(results)
    style = method_plot_style(results(i).name);
    semilogy(tEval, floor_log_values(results(i).orthLoss), ...
        'Color', style.Color, 'LineStyle', style.LineStyle, ...
        'LineWidth', style.LineWidth);
end
grid on;
set(gca, 'YScale', 'log');
xlabel('t');
ylabel('orthogonality loss');
legend({results.label}, 'Interpreter', 'none', 'Location', 'southoutside', 'NumColumns', 2);

if nargin >= 4 && ~isempty(outFile)
    hide_axes_toolbars(fig);
    exportgraphics(fig, outFile, 'Resolution', 300);
end
end


%% Source: .\matlab\common\projector_relative_error.m

function err = projector_relative_error(Uhat, Utrue)
%PROJECTOR_RELATIVE_ERROR Relative Frobenius projector error.

if any(~isfinite(Uhat(:))) || any(~isfinite(Utrue(:)))
    err = NaN;
    return;
end

p = size(Utrue, 2);
Qhat = orthonormalize(Uhat);
Qtrue = orthonormalize(Utrue);

% For two p-dimensional subspaces,
% ||Qhat Qhat' - Qtrue Qtrue'||_F = sqrt(2) ||(I-Qtrue Qtrue')Qhat||_F.
% This residual form avoids the severe cancellation in 2p-2||Qhat'Qtrue||_F^2
% when the two subspaces are nearly identical.
residual = Qhat - Qtrue*(Qtrue.'*Qhat);
err = sqrt(2)*norm(residual, 'fro') / sqrt(p);
end


%% Source: .\matlab\common\reconstruct_chart_point.m

function U = reconstruct_chart_point(chart, xrow)
%RECONSTRUCT_CHART_POINT Reconstruct a Stiefel representative from chart data.

p = chart.p;
q = chart.n - chart.p;
Xi = reshape(xrow(:), q, p);

if any(~isfinite(Xi(:)))
    U = NaN(chart.n, chart.p);
    return;
end

A = [eye(p); Xi];
try
    R = chol(eye(p) + Xi.'*Xi);
    V = A / R;
catch
    [V, ~] = thin_qr_positive(A);
end

U = apply_chart_inverse(chart, V);
end


%% Source: .\matlab\common\reconstruct_normal_point.m

function U = reconstruct_normal_point(coord, xrow)
%RECONSTRUCT_NORMAL_POINT Reconstruct a Stiefel representative from log data.

X = reshape(xrow(:), coord.n, coord.p);
if any(~isfinite(X(:)))
    U = NaN(coord.n, coord.p);
    return;
end
U = grassmann_exp(coord.Uref, X);
end


%% Source: .\matlab\common\recurrence_diagnostic.m

function cRec = recurrence_diagnostic(PhiAug, H)
%RECURRENCE_DIAGNOSTIC Simple a posteriori CV+A recurrence-growth diagnostic.

subdiag = abs(diag(H, -1));
if isempty(subdiag)
    cRec = max(vecnorm(PhiAug, 2, 2));
    return;
end

localAmp = zeros(numel(subdiag), 1);
for j = 1:numel(subdiag)
    localAmp(j) = (1 + sum(abs(H(1:j, j)))) / max(subdiag(j), realmin);
end

cRec = max(vecnorm(PhiAug, 2, 2)) * max(localAmp);
end


%% Source: .\matlab\common\repair_duplicate_rows.m

function selected = repair_duplicate_rows(A, selected)
%REPAIR_DUPLICATE_ROWS Replace duplicate selected rows by leverage candidates.

[n, p] = size(A);
if numel(unique(selected)) == p
    return;
end

used = false(n, 1);
fixed = zeros(1, p);
for j = 1:p
    idx = selected(j);
    if idx >= 1 && idx <= n && ~used(idx)
        fixed(j) = idx;
        used(idx) = true;
    end
end

scores = sum(abs(A), 2);
[~, order] = sort(scores, 'descend');
for j = 1:p
    if fixed(j) ~= 0
        continue;
    end
    replacement = order(find(~used(order), 1, 'first'));
    fixed(j) = replacement;
    used(replacement) = true;
end

selected = fixed;
end


%% Source: .\matlab\common\results_curve_table.m

function T = results_curve_table(results)
%RESULTS_CURVE_TABLE Long table of pointwise errors and orthogonality loss.

t = [];
method = {};
label = {};
coordType = {};
algebra = {};
err = [];
orth = [];

for i = 1:numel(results)
    nEval = numel(results(i).tEval);
    t = [t; results(i).tEval(:)]; %#ok<AGROW>
    method = [method; repmat({results(i).name}, nEval, 1)]; %#ok<AGROW>
    label = [label; repmat({results(i).label}, nEval, 1)]; %#ok<AGROW>
    coordType = [coordType; repmat({results(i).coordType}, nEval, 1)]; %#ok<AGROW>
    algebra = [algebra; repmat({results(i).algebra}, nEval, 1)]; %#ok<AGROW>
    err = [err; results(i).error(:)]; %#ok<AGROW>
    orth = [orth; results(i).orthLoss(:)]; %#ok<AGROW>
end

T = table(t, method, label, coordType, algebra, err, orth, ...
    'VariableNames', {'t', 'method', 'label', 'coordType', 'algebra', ...
    'relativeProjectorError', 'orthogonalityLoss'});
end


%% Source: .\matlab\common\results_summary_table.m

function T = results_summary_table(results)
%RESULTS_SUMMARY_TABLE Convert method summaries to a MATLAB table.

n = numel(results);
S = repmat(results(1).summary, n, 1);
for i = 1:n
    S(i) = results(i).summary;
end
T = struct2table(S);
end


%% Source: .\matlab\common\run_hermite_comparison.m

function out = run_hermite_comparison(tTrain, trainSample, tEval, trueSample, centerIndex, methods)
%RUN_HERMITE_COMPARISON Evaluate all selected Hermite methods.

if nargin < 6 || isempty(methods)
    methods = make_method_list('full');
end

results = [];
for i = 1:numel(methods)
    fprintf('  %-38s', methods(i).label);
    tic;
    result = evaluate_hermite_method(methods(i), tTrain, ...
        trainSample.U, trainSample.Delta, tEval, trueSample.U, centerIndex);
    if i == 1
        results = result;
    else
        results(i) = result; %#ok<AGROW>
    end
    fprintf(' max err = %.3e, time = %.2fs\n', ...
        results(i).summary.MaxError, toc);
end

out = struct();
out.results = results;
out.summaryTable = results_summary_table(results);
out.curveTable = results_curve_table(results);
end


%% Source: .\matlab\common\save_comparison_outputs.m

function save_comparison_outputs(cfg, prefix, comparison)
%SAVE_COMPARISON_OUTPUTS Write summary, pointwise curves, and MAT data.

if cfg.saveData
    writetable(comparison.summaryTable, ...
        fullfile(cfg.outputTableDir, [prefix, '_summary.csv']));
    writetable(comparison.curveTable, ...
        fullfile(cfg.outputTableDir, [prefix, '_curves.csv']));
    save(fullfile(cfg.outputDataDir, [prefix, '.mat']), 'comparison');
end
end


%% Source: .\matlab\common\select_results.m

function selected = select_results(results, methodNames)
%SELECT_RESULTS Keep results whose names occur in methodNames.

if nargin < 2 || isempty(methodNames)
    selected = results;
    return;
end

keep = false(numel(results), 1);
for i = 1:numel(results)
    keep(i) = any(strcmp(results(i).name, methodNames));
end
selected = results(keep);
end


%% Source: .\matlab\common\strip_large_coord.m

function light = strip_large_coord(coord)
%STRIP_LARGE_COORD Keep diagnostics but omit large coordinate arrays.

light = coord;
if isfield(light, 'X')
    light = rmfield(light, 'X');
end
if isfield(light, 'Xdot')
    light = rmfield(light, 'Xdot');
end
if isfield(light, 'chart') && isfield(light.chart, 'Q')
    light.chart.Q = [];
end
if isfield(light, 'Uref')
    light.Uref = [];
end
end


%% Source: .\matlab\common\summarize_method_result.m

function row = summarize_method_result(result, n, p, m)
%SUMMARIZE_METHOD_RESULT One-row diagnostics for a method result.

row = struct();
row.Method = result.name;
row.Label = result.label;
row.CoordType = result.coordType;
row.Algebra = result.algebra;
row.n = n;
row.p = p;
row.Nodes = m;
row.Degree = result.alg.degree;
row.MaxError = max(result.error, [], 'omitnan');
row.MedianError = median(result.error, 'omitnan');
row.MaxOrthogonalityLoss = max(result.orthLoss, [], 'omitnan');
row.ConfluentVandermondeCondition = result.alg.condSystem;
row.ArnoldiEtaQ = result.alg.etaQ;
row.ArnoldiCondQ = result.alg.condQ;
row.MinSubdiagonal = result.alg.minSubdiag;
row.BasisGrowth = result.alg.basisGrowth;
row.RecurrenceDiagnostic = result.alg.recurrenceDiagnostic;
row.CoefficientNormF = result.alg.coeffNormF;
row.Failure = result.failure;

if isfield(result.coord, 'kappaF')
    row.MaxChartKappaF = max(result.coord.kappaF, [], 'omitnan');
    row.MaxChartKappa2 = max(result.coord.kappa2, [], 'omitnan');
    row.MinChartSigma = min(result.coord.sigmaMin, [], 'omitnan');
else
    row.MaxChartKappaF = NaN;
    row.MaxChartKappa2 = NaN;
    row.MinChartSigma = NaN;
end

if isfield(result.coord, 'coordNorm')
    row.MaxCoordinateNormF = max(result.coord.coordNorm, [], 'omitnan');
else
    row.MaxCoordinateNormF = NaN;
end
end


%% Source: .\matlab\common\synthetic_curve_raw.m

function [Y, Ydot] = synthetic_curve_raw(curve, t)
%SYNTHETIC_CURVE_RAW Evaluate Y(t) and dY/dt for the synthetic curve.

Y = curve.Y0 + sin(3*t)*curve.Y1 + cos(3*t)*curve.Y2 + exp(t)*curve.Y3;
Ydot = 3*cos(3*t)*curve.Y1 - 3*sin(3*t)*curve.Y2 + exp(t)*curve.Y3;
end


%% Source: .\matlab\common\synthetic_curve_sample.m

function sample = synthetic_curve_sample(curve, t, opts)
%SYNTHETIC_CURVE_SAMPLE Stiefel representatives and horizontal derivatives.

if nargin < 3
    opts = struct();
end
if ~isfield(opts, 'noiseLevel')
    opts.noiseLevel = 0;
end
if ~isfield(opts, 'noiseSeed')
    opts.noiseSeed = 2718;
end

t = t(:);
m = numel(t);
U = cell(m, 1);
Delta = cell(m, 1);
Yraw = cell(m, 1);
Yused = cell(m, 1);
Ydot = cell(m, 1);

if opts.noiseLevel > 0
    state = rng;
    rng(opts.noiseSeed, 'twister');
end

for i = 1:m
    [Yi, YiDot] = synthetic_curve_raw(curve, t(i));
    if opts.noiseLevel > 0
        E = 2*rand(size(Yi)) - 1;
        E = E / max(norm(E, 'fro'), realmin);
        YiUsed = Yi + opts.noiseLevel*E;
    else
        YiUsed = Yi;
    end

    [Ui, Di] = horizontal_lift_from_full_rank(YiUsed, YiDot);
    U{i} = Ui;
    Delta{i} = Di;
    Yraw{i} = Yi;
    Yused{i} = YiUsed;
    Ydot{i} = YiDot;
end

if opts.noiseLevel > 0
    rng(state);
end

sample = struct();
sample.t = t;
sample.U = U;
sample.Delta = Delta;
sample.Yraw = Yraw;
sample.Yused = Yused;
sample.Ydot = Ydot;
sample.n = curve.n;
sample.p = curve.p;
sample.noiseLevel = opts.noiseLevel;
end


%% Source: .\matlab\common\synthetic_curve_setup.m

function curve = synthetic_curve_setup(n, p, seed)
%SYNTHETIC_CURVE_SETUP Random transcendental curve from the manuscript.

if nargin < 3
    seed = 1137;
end

state = rng;
rng(seed, 'twister');

curve = struct();
curve.n = n;
curve.p = p;
curve.seed = seed;
curve.Y0 = rand(n, p) - 0.5;
curve.Y1 = rand(n, p) - 0.5;
curve.Y2 = rand(n, p) - 0.5;
curve.Y3 = rand(n, p) - 0.5;

rng(state);
end


%% Source: .\matlab\common\thin_qr_positive.m

function [Q, R] = thin_qr_positive(Y)
%THIN_QR_POSITIVE Thin QR with positive diagonal in R when possible.

[Q, R] = qr(Y, 0);
d = sign(diag(R));
d(d == 0) = 1;
D = diag(d);
Q = Q*D;
R = D*R;
end


%% Source: .\matlab\common\va_arnoldi.m

function arn = va_arnoldi(t, degree)
%VA_ARNOLDI Scaled Arnoldi basis for Lagrange interpolation.

t = t(:);
m = numel(t);
if degree > m - 1
    error('V+A degree %d exceeds the Lagrange data dimension %d.', degree, m - 1);
end

Q = zeros(m, degree + 1);
H = zeros(degree + 1, degree);
Q(:, 1) = ones(m, 1);

for j = 1:degree
    v = t.*Q(:, j);
    for pass = 1:2
        for i = 1:j
            hij = (Q(:, i).'*v) / m;
            H(i, j) = H(i, j) + hij;
            v = v - Q(:, i)*hij;
        end
    end
    hnext = norm(v) / sqrt(m);
    H(j + 1, j) = hnext;
    if hnext <= 100*eps
        warning('va_arnoldi:Breakdown', ...
            'Near breakdown at step %d; remaining basis columns are set to NaN.', j);
        Q(:, j+1:end) = NaN;
        break;
    end
    Q(:, j + 1) = v / hnext;
end

arn = struct();
arn.Q = Q;
arn.H = H;
arn.degree = degree;
arn.m = m;
arn.etaQ = norm((Q.'*Q)/m - eye(degree + 1), 2);
arn.condQ = cond(Q);
arn.minSubdiag = min(abs(diag(H, -1)));
end


%% Source: .\matlab\common\vandermonde_values.m

function V = vandermonde_values(t, degree)
%VANDERMONDE_VALUES Monomial Vandermonde with columns 1,t,...,t^degree.

t = t(:);
V = ones(numel(t), degree + 1);
for j = 2:degree+1
    V(:, j) = V(:, j-1).*t;
end
end


%% Source: .\matlab\common\write_paper_figure_if_requested.m

function write_paper_figure_if_requested(cfg, localFigurePath, paperFigureName)
%WRITE_PAPER_FIGURE_IF_REQUESTED Copy generated figure to project root.

if cfg.writePaperFigures
    copyfile(localFigurePath, fullfile(cfg.projectRoot, paperFigureName));
end
end


%% Source: .\matlab\common\zero_tangent_cells.m

function Dcell = zero_tangent_cells(Ucell)
%ZERO_TANGENT_CELLS Create zero tangent cells matching Ucell sizes.

Dcell = cell(size(Ucell));
for i = 1:numel(Ucell)
    Dcell{i} = zeros(size(Ucell{i}));
end
end

