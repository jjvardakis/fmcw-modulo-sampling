%% EXPERIMENT 4 - Regime map: when does Hilbert-PoF beat plain unwrapping (Itoh)?
%  Noiseless, ideal IQ.  Random scenes: n targets on random bins 1..kmax, complex Gaussian amplitudes,
%  scaled so the peak of |rI|,|rQ| is P*lam.  Success = every fold recovered (max error < 1e-6).
%  Sufficient condition for Itoh:  max_n |r(n+1) - r(n)| < lam          (sample-to-sample change, per channel)
%  Ideal-model fold-count reference:   |M_C| < N/3                            (number of folds, both channels)
%  The empirical comparison tests where each method reconstructs successfully.
%  Peak levels are 1.25, 1.5 and 2 lam. Random scenes may lie on either
%  side of the step-size and fold-count reference conditions.
close all; setup_scene; rng(5);
quickRun = exist('v2_quick','var') && v2_quick;
T     = 100;                                  % trials per cell
ntar  = [1 2 4 8 16];                        % 1 x 5 number of targets
kmax  = [5 10 20 30 45 60 90 120];           % 1 x 8 highest bin a target may occupy (<= N/2-1)
peaks = [1.25 1.5 2];                        % 1 x 3 peak amplitude in units of lam
if quickRun
    T = 2;
    ntar = [1 4];
    kmax = [10 60];
    peaks = [1.25 2];
end

okI = zeros(numel(ntar), numel(kmax), numel(peaks));   % 5 x 8 x 3 Itoh success rate
okH = zeros(numel(ntar), numel(kmax), numel(peaks));   % 5 x 8 x 3 HPoF success rate
nMC = zeros(numel(ntar), numel(kmax), numel(peaks));   % 5 x 8 x 3 mean number of folds |M_C|
dAvg = zeros(numel(ntar), numel(kmax), numel(peaks));  % 5 x 8 x 3 mean of max|r(n+1)-r(n)| / lam
ex   = [];                                   % one example scene where Itoh fails (filled in the loop)
pts = zeros(0, 4);                           % one row per trial: [max|dr|/lam, |M_C|/(N/3), itoh ok, hpof ok]
for ip = 1:numel(peaks)
    for it = 1:numel(ntar)
        for ik = 1:numel(kmax)
            for t = 1:T
                k  = randperm(kmax(ik), min(ntar(it), kmax(ik)));          % 1 x n distinct bins in 1..kmax
                a  = randn(1, numel(k)) + 1i*randn(1, numel(k));           % 1 x n complex Gaussian amplitudes
                r  = fmcw_beat(k*dR, a, fc, B, Tc, N);                     % N x 1, targets on the range grid
                r  = r*peaks(ip)*lam/max(abs([real(r); imag(r)]));         % scale peak to P*lam
                zI = modlam(real(r), lam);  zQ = modlam(imag(r), lam);     % N x 1 each

                uok = max(abs(itoh(zI) + 1i*itoh(zQ) - r)) < 1e-6;          % Itoh exact?
                [hI, hQ] = hpof(zI, zQ, lam);
                hok = max(abs(hI + 1i*hQ - r)) < 1e-6;                      % HPoF exact?

                epsI = round((real(r) - zI)/(2*lam));  epsQ = round((imag(r) - zQ)/(2*lam));   % N x 1 integer residues
                MC   = numel(union(find(circshift(epsI,-1) ~= epsI), find(circshift(epsQ,-1) ~= epsQ)));
                dmax = max(abs([diff(real(r)); diff(imag(r))]));           % largest sample-to-sample change [V]

                okI(it, ik, ip) = okI(it, ik, ip) + uok/T;
                okH(it, ik, ip) = okH(it, ik, ip) + hok/T;
                nMC(it, ik, ip) = nMC(it, ik, ip) + MC/T;
                dAvg(it, ik, ip) = dAvg(it, ik, ip) + dmax/lam/T;
                if ~uok && isempty(ex) && ntar(it) == 2 && peaks(ip) == 1.5 && kmax(ik) == 30   % keep the first failing scene
                    ex = struct('r', r, 'zI', zI, 'itI', itoh(zI), 'dI', diff(real(r)));          % N x 1 each, (N-1) x 1
                end
                pts(end+1, :)   = [dmax/lam, MC/(N/3), uok, hok];           %#ok<SAGROW>  (#trials) x 4
            end
        end
    end
end

% ---- figure 1: ONE regime map per peak level (replaces the difference / folds / max-step columns) ----
%  Cell colour = bilinear mix of four corner colours, weighted by the two success rates pI (Itoh), pH (HPoF):
%     C = (1-pI)(1-pH)*Cfail + pI(1-pH)*Citoh + (1-pI)pH*Chpof + pI*pH*Cboth
%  These are blends of marginal success rates, not measured joint-outcome probabilities.
%  Cell text = both success rates in %.  Lines on cell edges = the two theoretical limits:
%     solid  : mean max|r(n+1)-r(n)| = lambda   (Itoh needs every step < lambda)
%     dotted : mean |M_C| = N/3                  (HPoF guarantee)
Cboth = [0.94 0.94 0.94];                    % 1 x 3 RGB  both succeed (light grey)
Citoh = [0.27 0.47 0.80];                    % 1 x 3 RGB  only Itoh    (blue)
Chpof = [0.85 0.33 0.24];                    % 1 x 3 RGB  only HPoF    (red)
Cfail = [0.22 0.22 0.24];                    % 1 x 3 RGB  both fail    (charcoal)
nt = numel(ntar);  nk = numel(kmax);
figure('Name', 'Exp4: regime map', 'Position', [30 30 1500 430]);
tiledlayout(1, numel(peaks), 'TileSpacing', 'compact', 'Padding', 'compact');
for ip = 1:numel(peaks)
    pI = okI(:,:,ip);  pH = okH(:,:,ip);     % nt x nk success rates (rows = targets, cols = kmax)
    img = zeros(nt, nk, 3);                  % nt x nk x 3 RGB image, one pixel per cell
    for c = 1:3
        img(:,:,c) = (1-pI).*(1-pH)*Cfail(c) + pI.*(1-pH)*Citoh(c) + (1-pI).*pH*Chpof(c) + pI.*pH*Cboth(c);
    end
    nexttile; image(1:nk, 1:nt, img); hold on              % image() puts row 1 (one target) at the top

    % thin white grid between cells
    for x = 1.5:1:nk-0.5, plot([x x], [0.5 nt+0.5], 'w-', 'LineWidth', 1); end
    for y = 1.5:1:nt-0.5, plot([0.5 nk+0.5], [y y], 'w-', 'LineWidth', 1); end

    % theoretical limits as staircases on the cell edges (black underlay + white line = visible on any colour)
    crit = {dAvg(:,:,ip) < 1, '-';  nMC(:,:,ip) < N/3, ':'};   % 2 x 2 cell: {nt x nk logical, line style}
    for ic = 1:2
        ok = crit{ic,1};
        for it = 1:nt
            for ik = 1:nk
                if ik < nk && ok(it,ik) ~= ok(it,ik+1)        % condition flips between this cell and the right one
                    xs = [ik ik] + 0.5;  ys = it + [-0.5 0.5];
                    plot(xs, ys, 'k-', 'LineWidth', 4);  plot(xs, ys, 'w', 'LineStyle', crit{ic,2}, 'LineWidth', 2);
                end
                if it < nt && ok(it,ik) ~= ok(it+1,ik)        % condition flips between this cell and the one below
                    xs = ik + [-0.5 0.5];  ys = [it it] + 0.5;
                    plot(xs, ys, 'k-', 'LineWidth', 4);  plot(xs, ys, 'w', 'LineStyle', crit{ic,2}, 'LineWidth', 2);
                end
            end
        end
    end

    % both success rates in every cell; white text on dark cells
    for it = 1:nt
        for ik = 1:nk
            rgb = squeeze(img(it, ik, :)).';                    % 1 x 3 colour of this cell
            tc  = 'k';  if rgb*[0.299; 0.587; 0.114] < 0.55, tc = 'w'; end
            text(ik, it, {sprintf('I %.0f', 100*pI(it,ik)), sprintf('H %.0f', 100*pH(it,ik))}, ...
                 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', tc);
        end
    end
    set(gca, 'XTick', 1:nk, 'XTickLabel', kmax, 'YTick', 1:nt, 'YTickLabel', ntar, 'TickLength', [0 0]);
    xlabel('highest target bin k_{max}');  title(sprintf('peak %g\\lambda', peaks(ip)));
    if ip == 1, ylabel('number of targets'); end
end
% one shared legend under all panels, built from invisible dummy objects (no colourbars needed)
h = [patch(NaN, NaN, Cboth), patch(NaN, NaN, Citoh), patch(NaN, NaN, Chpof), patch(NaN, NaN, Cfail), ...
     plot(NaN, NaN, 'k-', 'LineWidth', 2), plot(NaN, NaN, 'k:', 'LineWidth', 2)];   % 1 x 6 graphics handles
lgd = legend(h, {'both succeed', 'only Itoh', 'only HPoF', 'both fail', ...
                 'mean max|\Delta r| = \lambda (step reference)', 'mean |M_C| = N/3 (fold-count reference)'}, ...
             'Orientation', 'horizontal', 'NumColumns', 6, 'Box', 'off');
lgd.Layout.Tile = 'south';

% ---- printout: success per quadrant of figure 2 and per fold-ratio bin ----------------------
okU   = pts(:,3) == 1;  okP = pts(:,4) == 1;                              % (#trials) x 1 logical each
edges = [0 0.25 0.5 0.75 1 1.5 3 Inf];                                  % 1 x 8 bin edges of |M_C|/(N/3)
bin   = discretize(pts(:,2), edges);                                    % (#trials) x 1 bin index
rateU = accumarray(bin, okU, [numel(edges)-1 1], @mean);                % 7 x 1 Itoh success per bin
rateP = accumarray(bin, okP, [numel(edges)-1 1], @mean);                % 7 x 1 HPoF success per bin
cnt   = accumarray(bin, 1,   [numel(edges)-1 1]);                       % 7 x 1 trials per bin
q = {pts(:,1) < 1 & pts(:,2) < 1, pts(:,1) >= 1 & pts(:,2) < 1, pts(:,1) < 1 & pts(:,2) >= 1, pts(:,1) >= 1 & pts(:,2) >= 1};
qn = {'slow,  few folds ', 'fast,  few folds ', 'slow,  many folds', 'fast,  many folds'};
for i = 1:4
    fprintf('%s: %4d trials, Itoh %.2f, HPoF %.2f\n', qn{i}, nnz(q{i}), mean(okU(q{i})), mean(okP(q{i})));
end
fprintf('overall: Itoh %.3f, HPoF %.3f\n', mean(okU), mean(okP));
for i = 1:numel(cnt)
    fprintf('|M_C|/(N/3) in [%.2f,%.2f): %4d trials, Itoh %.2f, HPoF %.2f\n', edges(i), edges(i+1), cnt(i), rateU(i), rateP(i));
end