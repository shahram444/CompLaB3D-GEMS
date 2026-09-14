function exportSurrogateHeader(outFile, name, W, B, xoffset, xgain, xymin, ...
                               yoffset, ygain, yymin, trainMin, trainMax, ...
                               ylo, yhi, logOutput, meta, stats, microbeId, sourceCsv)
%EXPORTSURROGATEHEADER  Write a trained network as a CompLB3D C++ header.
%
%   This is the MATLAB twin of the emitter inside trainSurrogate.py.  The two
%   are kept BYTE-IDENTICAL on purpose: given the same weights they produce the
%   same file, so a MATLAB-trained network and a Python-trained one are
%   interchangeable and can be diffed against each other.  If you edit one
%   emitter, edit the other, and re-run testExporterParity.m to confirm.
%
%   Arguments
%     outFile    path of the .hh file to write
%     W          cell array of weight matrices, layer 1 first
%     B          cell array of bias column vectors, same length as W
%     xoffset    1 x nIn  mapminmax offsets  (the per-input minimum)
%     xgain      1 x nIn  mapminmax gains    (2 / range)
%     xymin      scalar, normally -1
%     yoffset, ygain, yymin   the same three for the output
%     trainMin, trainMax      1 x nIn, the box the data actually covered
%     ylo, yhi   observed growth range, 1/h
%     logOutput  true if the fit was done on log10(growth)
%     meta       struct with fields .model, .objective, .inputs (cell of strings)
%     stats      struct with fields .rmse_test, .r2_test, .maxerr_test
%     microbeId  global microbe index, for the pasted snippet
%     sourceCsv  training set filename, recorded in the header
%
%   Part of the CompLaB program.  GNU Affero General Public License v3 or later.

    nin   = size(W{1}, 2);
    guard = ['SURROGATE_WEIGHTS_' upper(name) '_HH'];
    ns    = ['srgnet_' name];

    f = fopen(outFile, 'w');
    if f < 0
        error('exportSurrogateHeader:open', 'Cannot open %s for writing.', outFile);
    end
    % Line endings are written explicitly as \n, not by using 'wt', so the file
    % is identical on Windows and Linux.  The C++ side does not care, but the
    % byte-identity check against the Python emitter does.
    w = @(s) fprintf(f, '%s\n', s);

    [~, csvBase, csvExt] = fileparts(sourceCsv);

    w('/* ==========================================================================');
    w(' * CompLaB surrogate network -- GENERATED FILE, DO NOT EDIT BY HAND.');
    w(' *');
    w(sprintf(' * Produced by trainSurrogate.m from: %s', [csvBase csvExt]));
    if ~isempty(meta.model)
        w(sprintf(' * Metabolic model              : %s', meta.model));
    end
    if ~isempty(meta.objective)
        w(sprintf(' * FBA objective                : %s', trunc(meta.objective, 66)));
    end
    w(' *');
    hidden = {};
    for i = 1:numel(W)-1
        hidden{end+1} = sprintf('%d tansig', size(W{i}, 1)); %#ok<AGROW>
    end
    w(sprintf(' * Architecture : %d input(s) -> %s -> 1 linear output', nin, strjoin(hidden, ' -> ')));
    np = 0;
    for i = 1:numel(W)
        np = np + numel(W{i}) + numel(B{i});
    end
    w(sprintf(' * Parameters   : %d weights and biases', np));
    if logOutput
        w(' * Output       : fitted in log space; eval() undoes the log for you.');
    end
    w(' *');
    w(' * FIT QUALITY (on data the network never saw during training)');
    w(sprintf(' *   test RMSE   : %s 1/h', trim6(stats.rmse_test)));
    w(sprintf(' *   test R^2    : %.6f', stats.r2_test));
    w(sprintf(' *   worst error : %s 1/h', trim6(stats.maxerr_test)));
    w(' *');
    w(' * TRAINING RANGE -- THE NETWORK IS ONLY VALID INSIDE THIS BOX.');
    w(' * A neural network extrapolates badly and will return confident');
    w(' * nonsense outside it.  inTrainingRange() below tests for exactly this;');
    w(' * call it if a run may wander out of range.');
    for k = 1:nin
        if k <= numel(meta.inputs)
            label = trunc(meta.inputs{k}, 58);
        else
            label = sprintf('input %d', k-1);
        end
        w(sprintf(' *   %-58s', label));
        w(sprintf(' *     %s .. %s mmol/gDW/h', trim10(trainMin(k)), trim10(trainMax(k))));
    end
    w(sprintf(' *   growth rate: %s .. %s 1/h', trim10(ylo), trim10(yhi)));
    w(' *');
    w(' * HOW TO USE IT -- two lines in surrogateModel.hh:');
    w(' *');
    w(sprintf(' *     #include "surrogate_weights_%s.hh"', name));
    w(' *');
    w(sprintf(' *     else if (microbeId == %d) {', microbeId));
    w(sprintf(' *         bioR[microbeId] = %s::eval(Fin[microbeId]);', ns));
    w(' *     }');
    w(' *');
    w(' * Fin is in mmol/gDW/h and POSITIVE means consumption; eval() returns the');
    w(' * specific growth rate in 1/h.  Those are CompLaB''s conventions already,');
    w(' * so no conversion is needed on either side.');
    w(' * ======================================================================== */');
    w(sprintf('#ifndef %s', guard));
    w(sprintf('#define %s', guard));
    w('');
    w('#include <vector>');
    w('#include <cmath>');
    w('#include <cstddef>');
    w('');
    w(sprintf('namespace %s {', ns));
    w('');
    w(sprintf('static const std::size_t nInputs = %d;', nin));
    w('');
    w('/* mapminmax input scaling: xs = (x - x_offset) * x_gain + x_ymin */');
    w(sprintf('static const double x_offset[%d] = { %s };', nin, joinNums(xoffset)));
    w(sprintf('static const double x_gain  [%d] = { %s };', nin, joinNums(xgain)));
    w(sprintf('static const double x_ymin      = %s;', num(xymin)));
    w('');
    w('/* mapminmax output scaling, applied in reverse: y = (ys - y_ymin)/y_gain + y_offset */');
    w(sprintf('static const double y_offset = %s;', num(yoffset)));
    w(sprintf('static const double y_gain   = %s;', num(ygain)));
    w(sprintf('static const double y_ymin   = %s;', num(yymin)));
    w('');
    w('/* the box the network was trained on */');
    w(sprintf('static const double trainMin[%d] = { %s };', nin, joinNums(trainMin)));
    w(sprintf('static const double trainMax[%d] = { %s };', nin, joinNums(trainMax)));
    w('');
    if logOutput
        w('static const bool logOutput = true;');
    else
        w('static const bool logOutput = false;');
    end
    w('');

    for i = 1:numel(W)
        Wi = W{i};
        Bi = B{i};
        [r, c] = size(Wi);
        if i == 1
            wname = 'IW1_1';
        else
            wname = sprintf('LW%d_%d', i, i-1);
        end
        bname = sprintf('b%d', i);
        w(sprintf('/* layer %d: %d x %d */', i, r, c));
        w(sprintf('static const double %s[%d][%d] = {', wname, r, c));
        for rr = 1:r
            w(sprintf('    { %s },', joinNums(Wi(rr, :))));
        end
        w('};');
        w(sprintf('static const double %s[%d] = { %s };', bname, r, joinNums(Bi(:).')));
        w('');
    end

    w('/* MATLAB''s tansig.  Algebraically identical to std::tanh; this is the');
    w(' * form MATLAB''s genFunction emits, kept so the generated file can be');
    w(' * compared line by line against a MATLAB-exported one. */');
    w('static inline double tansig(double n) { return 2.0 / (1.0 + std::exp(-2.0 * n)) - 1.0; }');
    w('');
    w('/* True when every input lies inside the trained box.  Outside it the');
    w(' * return value of eval() is not meaningful. */');
    w('static inline bool inTrainingRange(const std::vector<double>& x)');
    w('{');
    w('    if (x.size() < nInputs) return false;');
    w('    for (std::size_t i = 0; i < nInputs; ++i) {');
    w('        if (x[i] < trainMin[i] || x[i] > trainMax[i]) return false;');
    w('    }');
    w('    return true;');
    w('}');
    w('');
    w('/* Specific growth rate, 1/h.  x holds the uptake flux estimates in');
    w(' * mmol/gDW/h, POSITIVE for consumption, in CompLaB.xml substrate order. */');
    w('static inline double eval(const std::vector<double>& x)');
    w('{');
    w('    if (x.size() < nInputs) return 0.0;');
    w('');
    w(sprintf('    double a0[%d];', nin));
    w('    for (std::size_t i = 0; i < nInputs; ++i) {');
    w('        a0[i] = (x[i] - x_offset[i]) * x_gain[i] + x_ymin;');
    w('    }');
    w('');
    prev = 'a0';
    for i = 1:numel(W)
        [r, c] = size(W{i});
        cur = sprintf('a%d', i);
        if i == 1
            wname = 'IW1_1';
        else
            wname = sprintf('LW%d_%d', i, i-1);
        end
        bname = sprintf('b%d', i);
        w(sprintf('    double %s[%d];', cur, r));
        w(sprintf('    for (int r = 0; r < %d; ++r) {', r));
        w(sprintf('        double s = %s[r];', bname));
        w(sprintf('        for (int c = 0; c < %d; ++c) { s += %s[r][c] * %s[c]; }', c, wname, prev));
        if i == numel(W)
            w(sprintf('        %s[r] = s;', cur));
        else
            w(sprintf('        %s[r] = tansig(s);', cur));
        end
        w('    }');
        w('');
        prev = cur;
    end
    w(sprintf('    double g = (%s[0] - y_ymin) / y_gain + y_offset;', prev));
    if logOutput
        w('    /* the fit was done on log10(growth); undo it */');
        w('    g = std::pow(10.0, g);');
    end
    w('');
    w('    /* A tiny positive prediction is fit noise, not growth.  The shipped');
    w('     * network used the same 1e-8 cut-off; keeping it means a converted');
    w('     * network reproduces the old behaviour exactly. */');
    w('    if (!(g > 1e-8)) { g = 0.0; }   /* also catches NaN */');
    w('    return g;');
    w('}');
    w('');
    w(sprintf('}  // namespace %s', ns));
    w('');
    w(sprintf('#endif  // %s', guard));

    fclose(f);
    fprintf('Wrote %s\n', outFile);
end


% --- helpers ---------------------------------------------------------------

function s = num(v)
    % 17 significant digits round-trips an IEEE double exactly.  Anything
    % shorter silently changes the network.
    s = sprintf('%.17g', v);
end

function s = joinNums(v)
    parts = cell(1, numel(v));
    for i = 1:numel(v)
        parts{i} = sprintf('%.17g', v(i));
    end
    s = strjoin(parts, ', ');
end

function s = trim6(v)
    s = sprintf('%.6g', v);
end

function s = trim10(v)
    s = sprintf('%.10g', v);
end

function s = trunc(str, n)
    if numel(str) > n
        s = str(1:n);
    else
        s = str;
    end
end
