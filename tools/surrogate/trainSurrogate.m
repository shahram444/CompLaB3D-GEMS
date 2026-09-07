function trainSurrogate(csvFile, varargin)
%TRAINSURROGATE  Train a CompLaB surrogate network and export it as C++.
%
%   This is the MATLAB path, using the Deep Learning Toolbox, and it
%   reconstructs the workflow that produced the network originally shipped in
%   surrogateModel.hh.  That network's variable names (IW1_1, LW2_1, tansig,
%   the mapminmax gains and offsets) are MATLAB genFunction output, but the
%   training script itself was never distributed with CompLaB, so the weights
%   could be evaluated and never reproduced.  This file closes that gap.
%
%   trainSurrogate.py does the same job with numpy/scipy and no licence, and
%   writes a byte-identical header.  Use whichever suits you.
%
%   USAGE
%       trainSurrogate('geobacter_training.csv', ...
%                      'Name', 'geobacter', ...
%                      'MicrobeId', 0, ...
%                      'Output', 'surrogate_weights_geobacter.hh')
%
%   OPTIONS  (name/value)
%       'Name'        identifier for the namespace and include guard  ['surrogate']
%       'MicrobeId'   global microbe index, as ordered in CompLaB.xml  [0]
%       'Layers'      hidden layer sizes                        [[10 10 10 10]]
%       'Restarts'    independent fits; best validation error wins        [5]
%       'LogOutput'   fit log10(growth) instead of growth             [false]
%       'Seed'        RNG seed, for reproducibility                       [0]
%       'Output'      output .hh path             ['surrogate_weights_<Name>.hh']
%
%   WHY THE DEFAULTS ARE WHAT THEY ARE
%       [10 10 10 10] tansig hidden layers with a linear output, mapminmax on
%       both ends, is exactly the shape of the shipped network.  Keeping it as
%       the default means a retrained network drops into an existing build
%       with no other change.
%
%       trainlm (Levenberg-Marquardt) is MATLAB's default for small nets and
%       is the right choice here: a few hundred parameters, a smooth response
%       surface, and a few thousand samples.  Early stopping on a 15%%
%       validation split is what stops it memorising the FBA sweep.
%
%   Part of the CompLaB program.  GNU Affero General Public License v3 or later.

    p = inputParser;
    p.addParameter('Name',      'surrogate', @ischar);
    p.addParameter('MicrobeId', 0,           @isnumeric);
    p.addParameter('Layers',    [10 10 10 10], @isnumeric);
    p.addParameter('Restarts',  5,           @isnumeric);
    p.addParameter('LogOutput', false,       @islogical);
    p.addParameter('Seed',      0,           @isnumeric);
    p.addParameter('Output',    '',          @ischar);
    p.parse(varargin{:});
    opt = p.Results;

    if isempty(opt.Output)
        opt.Output = sprintf('surrogate_weights_%s.hh', opt.Name);
    end

    if isempty(which('feedforwardnet'))
        error('trainSurrogate:noToolbox', ...
              ['The Deep Learning Toolbox is not on the path (feedforwardnet ' ...
               'not found).\nUse the Python path instead:\n' ...
               '    python3 trainSurrogate.py %s -o %s'], csvFile, opt.Output);
    end

    % ---- load ------------------------------------------------------------
    [X, y, meta, header] = loadTrainingCsv(csvFile);
    nin = size(X, 2);
    fprintf('Loaded %s: %d samples, %d input(s)\n', csvFile, size(X,1), nin);
    if ~isempty(header)
        fprintf('  columns: %s\n', strjoin(header, ', '));
    end

    trainMin = min(X, [], 1);
    trainMax = max(X, [], 1);
    ylo = min(y);
    yhi = max(y);
    fprintf('  input range : %s\n', strjoin(arrayfun(@(a,b) sprintf('[%.6g, %.6g]', a, b), ...
                                            trainMin, trainMax, 'UniformOutput', false), ', '));
    fprintf('  growth range: [%.6g, %.6g] 1/h\n', ylo, yhi);

    if yhi <= 0
        error('trainSurrogate:allZero', ...
              ['Every sample has zero growth; there is nothing to fit. ' ...
               'Re-run generateTrainingData.py over a range where the ' ...
               'organism actually grows.']);
    end

    yfit = y;
    if opt.LogOutput
        floorv = max(yhi * 1e-9, 1e-12);
        nclip  = sum(yfit < floorv);
        yfit   = log10(max(yfit, floorv));
        fprintf('  fitting log10(growth); %d zero/near-zero sample(s) floored at %.3g\n', ...
                nclip, floorv);
    end

    % MATLAB's neural nets take COLUMNS as samples.  Getting this transpose
    % wrong is the single most common way to produce a net that trains to a
    % flat line without any error message.
    Xt = X.';
    Yt = yfit(:).';

    % ---- train -----------------------------------------------------------
    best = [];
    bestVal = inf;
    for k = 1:opt.Restarts
        rng(opt.Seed + 1000*k);

        net = feedforwardnet(opt.Layers, 'trainlm');

        % Force the processing chain to mapminmax ONLY.  By default MATLAB
        % also inserts 'removeconstantrows', which silently deletes an input
        % column if the sweep happened to hold it constant -- and then the
        % exported weight matrix has fewer columns than CompLaB will pass in,
        % and the C++ reads past the end of the row.
        net.inputs{1}.processFcns  = {'mapminmax'};
        net.outputs{end}.processFcns = {'mapminmax'};

        for i = 1:numel(opt.Layers)
            net.layers{i}.transferFcn = 'tansig';
        end
        net.layers{end}.transferFcn = 'purelin';   % linear output

        net.divideFcn   = 'dividerand';
        net.divideParam.trainRatio = 0.70;
        net.divideParam.valRatio   = 0.15;
        net.divideParam.testRatio  = 0.15;
        net.trainParam.showWindow  = false;
        net.trainParam.epochs      = 1000;
        net.trainParam.max_fail    = 20;           % early-stopping patience

        [net, tr] = train(net, Xt, Yt);
        vperf = tr.best_vperf;
        marker = '';
        if vperf < bestVal
            marker = '   <- best so far';
        end
        fprintf('  restart %d/%d: train MSE %.4e, validation MSE %.4e%s\n', ...
                k, opt.Restarts, tr.best_perf, vperf, marker);
        if vperf < bestVal
            bestVal = vperf;
            best = struct('net', net, 'tr', tr);
        end
    end
    net = best.net;
    tr  = best.tr;

    % ---- pull the weights and the mapminmax settings out -----------------
    nlayers = numel(opt.Layers) + 1;
    W = cell(1, nlayers);
    B = cell(1, nlayers);
    W{1} = net.IW{1,1};
    B{1} = net.b{1};
    for i = 2:nlayers
        W{i} = net.LW{i, i-1};
        B{i} = net.b{i};
    end

    ps = net.inputs{1}.processSettings{1};
    xoffset = ps.xoffset(:).';
    xgain   = ps.gain(:).';
    xymin   = ps.ymin;

    qs = net.outputs{end}.processSettings{1};
    yoffset = qs.xoffset;
    ygain   = qs.gain;
    yymin   = qs.ymin;

    % ---- fit quality on the held-out test split --------------------------
    pred = @(Xraw) unlog(mapminmaxReverse(sim(net, Xraw.'), yoffset, ygain, yymin), opt.LogOutput);
    stats = struct();
    splits = {'train', tr.trainInd; 'val', tr.valInd; 'test', tr.testInd};
    fprintf('\nFit quality\n');
    for s = 1:size(splits, 1)
        idx = splits{s, 2};
        pv  = pred(X(idx, :));
        pv(pv <= 1e-8) = 0;
        e   = pv(:) - y(idx);
        ss  = sum((y(idx) - mean(y(idx))).^2);
        r2  = 1 - (e.'*e)/ss;
        stats.(['rmse_'   splits{s,1}]) = sqrt(mean(e.^2));
        stats.(['maxerr_' splits{s,1}]) = max(abs(e));
        stats.(['r2_'     splits{s,1}]) = r2;
        fprintf('  %-5s  RMSE %.6g 1/h   max err %.6g 1/h   R^2 %.6f\n', ...
                splits{s,1}, sqrt(mean(e.^2)), max(abs(e)), r2);
    end

    if stats.rmse_test / max(yhi, realmin) > 0.05
        warning('trainSurrogate:poorFit', ...
                ['Test RMSE is %.1f%% of the largest growth rate in the ' ...
                 'training set. That is a poor fit. Try more samples, more ' ...
                 'restarts, ''LogOutput'' if growth spans orders of ' ...
                 'magnitude, or a narrower input range.'], ...
                100 * stats.rmse_test / yhi);
    end

    % ---- export ----------------------------------------------------------
    exportSurrogateHeader(opt.Output, opt.Name, W, B, xoffset, xgain, xymin, ...
                          yoffset, ygain, yymin, trainMin, trainMax, ...
                          ylo, yhi, opt.LogOutput, meta, stats, ...
                          opt.MicrobeId, csvFile);

    % ---- reference table, so verifyExport.py can check the C++ -----------
    refFile = [regexprep(opt.Output, '\.hh$', '') '_reference.csv'];
    rng(12345);
    Xr = trainMin + rand(512, nin) .* (trainMax - trainMin);
    pr = pred(Xr);
    pr(pr <= 1e-8) = 0;
    fid = fopen(refFile, 'w');
    fprintf(fid, '# reference predictions for %s -- feed these to the compiled\n', opt.Output);
    fprintf(fid, '# header and the outputs must agree to ~1e-12 relative.\n');
    for i = 1:size(Xr, 1)
        fprintf(fid, '%s,%.17g\n', strjoin(arrayfun(@(v) sprintf('%.17g', v), ...
                Xr(i,:), 'UniformOutput', false), ','), pr(i));
    end
    fclose(fid);
    fprintf('Wrote %s (512 checkpoints for verifyExport.py)\n', refFile);

    fprintf('\nAdd these two things to surrogateModel.hh:\n\n');
    fprintf('    #include "%s"\n\n', opt.Output);
    fprintf('    else if (microbeId == %d) {\n', opt.MicrobeId);
    fprintf('        bioR[microbeId] = srgnet_%s::eval(Fin[microbeId]);\n', opt.Name);
    fprintf('    }\n\n');
    fprintf('and set <reaction_type>surrogate</reaction_type> plus\n');
    fprintf('<enable_surrogate>true</enable_surrogate> in CompLaB.xml.\n');
end


% ---------------------------------------------------------------------------
function y = mapminmaxReverse(ys, xoffset, gain, ymin)
    y = (ys - ymin) ./ gain + xoffset;
end

function y = unlog(y, doLog)
    if doLog
        y = 10.^y;
    end
end


function [X, y, meta, header] = loadTrainingCsv(path)
%LOADTRAININGCSV  Read generateTrainingData.py's CSV, metadata included.
    fid = fopen(path, 'r');
    if fid < 0
        error('trainSurrogate:open', 'Cannot open %s', path);
    end
    meta = struct('model', '', 'objective', '', 'inputs', {{}});
    header = {};
    rows = [];
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        s = strtrim(line);
        if isempty(s), continue; end
        if s(1) == '#'
            body = strtrim(s(2:end));
            if strncmp(body, 'model:', 6)
                meta.model = strtrim(body(7:end));
            elseif strncmp(body, 'objective:', 10)
                meta.objective = strtrim(body(11:end));
            elseif strncmp(body, 'input ', 6)
                meta.inputs{end+1} = body;
            end
            continue;
        end
        parts = strsplit(s, ',');
        if isempty(header) && any(isletter(parts{1}))
            header = strtrim(parts);
            continue;
        end
        rows(end+1, :) = cellfun(@str2double, parts); %#ok<AGROW>
    end
    fclose(fid);
    if isempty(rows)
        error('trainSurrogate:empty', '%s contains no data rows.', path);
    end
    if size(rows, 2) < 2
        error('trainSurrogate:tooFewColumns', ...
              ['%s needs at least one input column and one output column; ' ...
               'found %d.'], path, size(rows, 2));
    end
    X = rows(:, 1:end-1);
    y = rows(:, end);
end
