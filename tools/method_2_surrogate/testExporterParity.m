% TESTEXPORTERPARITY  Prove the MATLAB and Python emitters agree byte for byte.
%
% There are two independent implementations of the C++ emitter -- one in
% trainSurrogate.py, one in exportSurrogateHeader.m -- because not everyone has
% a MATLAB licence and not everyone wants a Python toolchain.  Two
% implementations of the same format is exactly the situation where they drift
% apart silently, so this script keeps them honest:
%
%   1. read a header the PYTHON emitter produced,
%   2. pull the weights and scalings back out of it,
%   3. re-emit with the MATLAB emitter,
%   4. diff.
%
% The only line allowed to differ is the provenance line naming which script
% produced the file.  Anything else is a bug in one of the two emitters.
%
% This test needs no toolbox and runs in GNU Octave as well as MATLAB, so it
% can be run in CI.
%
%   Usage:   testExporterParity('surrogate_weights_ecolicore.hh')
%
% Part of the CompLaB program.  GNU Affero General Public License v3 or later.

function testExporterParity(pyHeader)
    if nargin < 1
        error('usage: testExporterParity(''surrogate_weights_<name>.hh'')');
    end
    txt = fileread(pyHeader);

    name = regexp(txt, 'namespace srgnet_(\w+)', 'tokens', 'once');
    name = name{1};

    nin = str2double(getTok(txt, 'nInputs = (\d+);'));

    xoffset  = getVec(txt, 'x_offset');
    xgain    = getVec(txt, 'x_gain  ');
    xymin    = str2double(getTok(txt, 'x_ymin\s*=\s*(\S+);'));
    yoffset  = str2double(getTok(txt, 'y_offset\s*=\s*(\S+);'));
    ygain    = str2double(getTok(txt, 'y_gain\s*=\s*(\S+);'));
    yymin    = str2double(getTok(txt, 'y_ymin\s*=\s*(\S+);'));
    trainMin = getVec(txt, 'trainMin');
    trainMax = getVec(txt, 'trainMax');
    logOutput = ~isempty(strfind(txt, 'logOutput = true'));

    % weights, layer by layer, until a layer name is missing
    W = {}; B = {};
    i = 1;
    while true
        if i == 1
            wname = 'IW1_1';
        else
            wname = sprintf('LW%d_%d', i, i-1);
        end
        M = getMat(txt, wname);
        if isempty(M), break; end
        W{i} = M;
        B{i} = getVec(txt, sprintf('b%d', i));
        i = i + 1;
    end
    fprintf('Recovered %d layers from %s: ', numel(W), pyHeader);
    for k = 1:numel(W)
        fprintf('%dx%d ', size(W{k},1), size(W{k},2));
    end
    fprintf('\n');
    assert(size(W{1},2) == nin, 'input width does not match nInputs');

    % metadata and fit statistics, read back out of the comment block
    meta = struct();
    meta.model     = strtrimSafe(getTok(txt, 'Metabolic model\s+: ([^\n]*)'));
    meta.objective = strtrimSafe(getTok(txt, 'FBA objective\s+: ([^\n]*)'));
    meta.inputs = {};
    lines = strsplit(txt, sprintf('\n'), 'CollapseDelimiters', false);
    for k = 1:numel(lines)
        t = regexp(lines{k}, '^\s*\*\s+(input \d+: .*?)\s*$', 'tokens', 'once');
        if ~isempty(t)
            meta.inputs{end+1} = t{1};
        end
    end

    stats = struct();
    stats.rmse_test   = str2double(getTok(txt, 'test RMSE\s+: (\S+) 1/h'));
    stats.r2_test     = str2double(getTok(txt, 'test R\^2\s+: (\S+)'));
    stats.maxerr_test = str2double(getTok(txt, 'worst error\s+: (\S+) 1/h'));

    rangeLine = regexp(txt, 'growth rate: (\S+) \.\. (\S+) 1/h', 'tokens', 'once');
    ylo = str2double(rangeLine{1});
    yhi = str2double(rangeLine{2});

    microbeId = str2double(getTok(txt, 'else if \(microbeId == (\d+)\)'));
    srcCsv    = strtrimSafe(getTok(txt, 'Produced by trainSurrogate\.\w+ from: ([^\n]*)'));

    out = [tempname() '.hh'];
    exportSurrogateHeader(out, name, W, B, xoffset, xgain, xymin, ...
                          yoffset, ygain, yymin, trainMin, trainMax, ...
                          ylo, yhi, logOutput, meta, stats, microbeId, srcCsv);

    a = strsplit(fileread(pyHeader), sprintf('\n'), 'CollapseDelimiters', false);
    b = strsplit(fileread(out),      sprintf('\n'), 'CollapseDelimiters', false);

    % normalise the one line that is SUPPOSED to differ
    a = regexprep(a, 'trainSurrogate\.py from:', 'trainSurrogate.X from:');
    b = regexprep(b, 'trainSurrogate\.m from:',  'trainSurrogate.X from:');

    nbad = 0;
    n = max(numel(a), numel(b));
    for k = 1:n
        la = ''; lb = '';
        if k <= numel(a), la = a{k}; end
        if k <= numel(b), lb = b{k}; end
        if ~strcmp(la, lb)
            nbad = nbad + 1;
            if nbad <= 10
                fprintf('line %d differs:\n  py: %s\n  m : %s\n', k, la, lb);
            end
        end
    end
    if nbad == 0
        fprintf('PASS: %d lines, byte-identical between the two emitters.\n', numel(a));
    else
        error('FAIL: %d line(s) differ between the Python and MATLAB emitters.', nbad);
    end
end


function s = getTok(txt, pat)
    t = regexp(txt, pat, 'tokens', 'once');
    if isempty(t)
        s = '';
    else
        s = t{1};
    end
end

function s = strtrimSafe(s)
    if isempty(s)
        s = '';
    else
        s = strtrim(s);
    end
end

function v = getVec(txt, name)
    pat = [regexptranslate('escape', name) '\s*(?:\[\d+\])?\s*=\s*\{([^}]*)\};'];
    t = regexp(txt, pat, 'tokens', 'once');
    if isempty(t)
        v = [];
        return;
    end
    v = str2double(strsplit(strtrim(t{1}), ','));
end

function M = getMat(txt, name)
    pat = [regexptranslate('escape', name) '\[\d+\]\[\d+\] = \{(.*?)\n\};'];
    t = regexp(txt, pat, 'tokens', 'once');
    if isempty(t)
        M = [];
        return;
    end
    rows = regexp(t{1}, '\{([^{}]*)\}', 'tokens');
    M = [];
    for k = 1:numel(rows)
        M(k, :) = str2double(strsplit(strtrim(rows{k}{1}), ',')); %#ok<AGROW>
    end
end
