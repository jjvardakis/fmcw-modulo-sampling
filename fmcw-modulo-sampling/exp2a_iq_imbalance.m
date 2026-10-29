%% EXPERIMENT 2a - Separate gain and phase imbalance sweeps
% Compare ordinary HPoF with HPoF using the known receiver mismatch.
% Scene parameters and the NMSE helper come from setup_scene.m.
% Success means reconstruction NMSE < -100 dB relative to the imbalanced input.
% It does not necessarily mean that the selected model order is correct.

close all;
setup_scene;

% This experiment uses the report's 8 m / 30 m scene.
R = [8 30];
amp = amp_strong * (R(1) ./ R).^2;
kbin = round(R / dR);
rng(2);

quickRun = exist('v2_quick', 'var') && v2_quick;

%% Sweep settings

T = 50;                         % Trials per point

gdB = 0:0.25:3;                 % Gain error [dB]
phdeg = 0:5:30;                 % Phase error [degrees]

success_threshold_dB = -100;

% Small execution check, consistent with the other experiments.
if quickRun
    T = 2;
    gdB = [0 0.5 1];
    phdeg = [0 5 20];
end

% Each row contains [gain increment g, phase error phi].
% hpof uses a Q-channel gain of 1 + g.
% Gain and phase are varied separately: gain first, then phase.
sweep = [10.^(gdB(:) / 20) - 1, zeros(numel(gdB), 1); ...
         zeros(numel(phdeg), 1), phdeg(:) * pi / 180];

nmse_naive = zeros(size(sweep, 1), T);
nmse_aware = zeros(size(sweep, 1), T);

%% Run the trials

for i = 1:size(sweep, 1)

    g = sweep(i, 1);
    phi = sweep(i, 2);

    for trial = 1:T

        % Generate the FMCW beat signal with random target phases.
        r = fmcw_beat(R, amp .* exp(1i * 2 * pi * rand(1, numel(R))), ...
            fc, B, Tc, N);

        % Apply receiver imbalance to Q; I is unchanged.
        xI = real(r);
        xQ = (1 + g) * (imag(r) * cos(phi) + real(r) * sin(phi));

        % The reconstruction reference is the input before modulo folding.
        ref = xI + 1i * xQ;

        zI = modlam(xI, lam);
        zQ = modlam(xQ, lam);

        % Ordinary HPoF assumes ideal I/Q.
        [aI, aQ] = hpof(zI, zQ, lam);

        % Known-mismatch HPoF receives the true gain and phase errors.
        [bI, bQ] = hpof(zI, zQ, lam, g, phi);

        nmse_naive(i, trial) = nmse(aI + 1i * aQ, ref);
        nmse_aware(i, trial) = nmse(bI + 1i * bQ, ref);

    end
end

%% Calculate reconstruction success rates

success_naive = 100 * mean(nmse_naive < success_threshold_dB, 2);
success_aware = 100 * mean(nmse_aware < success_threshold_dB, 2);

ig = 1:numel(gdB);
ip = numel(gdB) + (1:numel(phdeg));

%% Plot the two sweeps side by side

f = figure('Name', 'Exp2a: IQ imbalance success - two panels', ...
    'Color', 'w', 'Position', [50 50 1100 420]);

subplot(1, 2, 1);

plot(gdB, success_naive(ig), 'bo-', ...
     gdB, success_aware(ig), 'g^--', 'LineWidth', 1.3);

grid on;
xlabel('Gain error [dB] (phase error = 0)');
ylabel('Successful reconstructions [%]');
ylim([-3 103]);
yticks(0:20:100);

legend('Ordinary HPoF', 'HPoF with known mismatch', ...
    'Location', 'southwest');

title('Gain-error sweep');

subplot(1, 2, 2);

plot(phdeg, success_naive(ip), 'bo-', ...
     phdeg, success_aware(ip), 'g^--', 'LineWidth', 1.3);

grid on;
xlabel('Phase error [degrees] (gain error = 0 dB)');
ylabel('Successful reconstructions [%]');
ylim([-3 103]);
yticks(0:20:100);
title('Phase-error sweep');

sgtitle(sprintf( ...
    'Success: NMSE < %g dB; %d trials per point; N = %d; ranges %s m', ...
    success_threshold_dB, T, N, mat2str(R)));

%% Display the measured success rates
% run_experiments saves figures and numerical arrays for all experiments.
fprintf('Scene ranges: %s m\n', mat2str(R));
fprintf('Gain sweep: gain dB, ordinary success percent, known-mismatch success percent\n');
disp([gdB(:), success_naive(ig), success_aware(ig)]);
fprintf('Phase sweep: phase degrees, ordinary success percent, known-mismatch success percent\n');
disp([phdeg(:), success_naive(ip), success_aware(ip)]);
