% PA_poly_regrowth.m  --  Memoryless Polynomial PA Spectral Regrowth
% =============================================================================
% Purpose:
%   Visualize per-order spectral regrowth of a memoryless polynomial PA:
%       y = sum_{k=1}^{K} a_k * x * |x|^(k-1)
%   For each order k, compute and plot the individual basis term's spectrum.
%
%   (a) Odd orders: k = 1, 3, 5, 7, 9
%   (b) All orders: k = 1, 2, 3, 4, 5, 6, 7, 8, 9
%
% Key question:
%   Do even-order terms contribute to ACLR?
%   (Two-part answer: physical RF vs baseband simulation -- see section 10)
%
% Input signal:
%   Complex baseband, bandlimited random with flat 20 MHz rectangular PSD,
%   generated via random-phase IFFT of a rectangular mask.  Fs = 240 MHz.
%
% Coefficients:
%   a_k = 1 for all k (analysis focuses on shape/position, not PA fit).
% =============================================================================

clear; clc; close all;

%% 1. Parameters
Fs    = 240e6;        % sampling rate [Hz]
BW    = 20e6;         % signal bandwidth [Hz]
N     = 2^15;         % FFT length (32768), df ~ 7.32 kHz
K_max = 9;            % highest polynomial order
rng(0);               % reproducibility

df   = Fs/N;
f    = ((-N/2):(N/2-1))' * df;       % freq axis (Hz), centered at 0
fMHz = f/1e6;

fprintf('=== Memoryless Polynomial PA -- Per-Order Spectral Regrowth ===\n');
fprintf('Fs = %.1f MHz,  BW = %.1f MHz,  N = %d,  df = %.2f kHz\n', ...
        Fs/1e6, BW/1e6, N, df/1e3);
fprintf('Coefficients a_k = 1 (all orders), basis y_k = x .* |x|^(k-1)\n\n');

%% 2. Signal Generation -- bandlimited rect via random-phase IFFT
mask  = double(abs(f) <= BW/2);
nbin  = nnz(mask);
phase = exp(1j * 2*pi * rand(N,1));
X_in  = mask .* phase;                   % complex rectangular spectrum

x_raw = ifft(ifftshift(X_in));           % complex baseband signal
x     = x_raw / max(abs(x_raw));         % peak-normalize to 1

rms_x  = sqrt(mean(abs(x).^2));
papr_dB = 20*log10(1/rms_x);
fprintf('Signal: complex baseband, %d active bins\n', nbin);
fprintf('        peak = 1, RMS = %.4f, PAPR = %.2f dB\n\n', rms_x, papr_dB);

%% 3. Input Spectrum Verification + Figure 0
X_check = fftshift(fft(x)) / N;
Xdb     = 20*log10(abs(X_check) + 1e-20);
Xdb     = Xdb - max(Xdb);                    % peak -> 0 dB

% in-band ripple (exclude edge bins)
ib_idx    = abs(f) <= BW/2 - 3*df;
ripple_dB = max(Xdb(ib_idx)) - min(Xdb(ib_idx));

% out-of-band rejection
oob_idx = abs(f) >= BW/2 + 5*df;
oob_max = max(Xdb(oob_idx));

fprintf('Input in-band ripple  : %.3f dB\n', ripple_dB);
fprintf('Input OOB max level   : %.1f dB\n\n', oob_max);

figure('Position',[80 80 1100 420]);
subplot(1,2,1);
plot(fMHz, Xdb, 'b-', 'LineWidth', 1.4); hold on;
xline(-BW/2/1e6, 'k--', 'LineWidth', 1);
xline( BW/2/1e6, 'k--', 'LineWidth', 1);
xlabel('Frequency (MHz)');
ylabel('|X(f)|  (dB, peak-normalized)');
title('Input Spectrum -- 20 MHz Rectangle');
xlim([-120 120]); ylim([-120 5]); grid on;

subplot(1,2,2);
histogram(abs(x), 80, 'FaceColor', [0.2 0.5 1.0], 'EdgeColor','none');
xlabel('|x|'); ylabel('count');
title(sprintf('Input |x| distribution   (PAPR = %.2f dB)', papr_dB));
grid on;

sgtitle('Input Signal Verification -- Bandlimited Rectangular PSD');
saveas(gcf, 'PA_poly_input.png');
fprintf('Saved PA_poly_input.png\n');

%% 4. Per-Order Basis Computation
% Hann window with ENBW normalization (preserves bin power)
w = hann(N, 'periodic');
w = w / sqrt(mean(w.^2));

Pk_lin = zeros(N, K_max);    % |Y_k(f)|^2
yk_td  = zeros(N, K_max);    % time-domain y_k
for k = 1:K_max
    yk           = x .* abs(x).^(k-1);
    yk_td(:,k)   = yk;
    Yk           = fftshift(fft(yk .* w)) / N;
    Pk_lin(:,k)  = abs(Yk).^2;
end

%% 5. Numerical Verification
fprintf('\n--- Per-order verification ---\n');
fprintf('%5s %15s %15s %15s\n', ...
        'Order','Parseval rel','Meas BW (MHz)','Expected(MHz)');
for k = 1:K_max
    % Parseval check (no window)
    Yk_nowin = fft(yk_td(:,k));
    lhs = sum(abs(yk_td(:,k)).^2);
    rhs = sum(abs(Yk_nowin).^2)/N;
    p_err = abs(lhs - rhs)/lhs;

    % -40 dB bandwidth
    Pk_norm   = Pk_lin(:,k) / max(Pk_lin(:,k));
    idx_above = find(10*log10(Pk_norm + 1e-30) > -40);
    if isempty(idx_above)
        bw_meas_MHz = 0;
    else
        bw_meas_MHz = (fMHz(max(idx_above)) - fMHz(min(idx_above)));
    end

    fprintf('%5d %15.2e %15.1f %15.1f\n', k, p_err, bw_meas_MHz, k*BW/1e6);
end

% Alias safety: k=9 energy near Nyquist should be below peak
k9_peak     = max(Pk_lin(:,K_max));
nyquist_idx = abs(fMHz) > 100;
if any(nyquist_idx)
    k9_nyq   = max(Pk_lin(nyquist_idx, K_max));
    alias_db = 10*log10(k9_nyq/k9_peak);
else
    alias_db = -Inf;
end
fprintf('k=9 energy at |f|>100 MHz : %.1f dBc (should be <0)\n', alias_db);

%% 6. Figure 1 -- Odd-Order Overlay (normalized to own peak)
odd_orders = [1 3 5 7 9];
colors_odd = [0.00 0.45 0.74;
              0.85 0.33 0.10;
              0.93 0.69 0.13;
              0.49 0.18 0.56;
              0.47 0.67 0.19];

figure('Position',[80 80 1100 600]);
hold on;
for i = 1:length(odd_orders)
    k   = odd_orders(i);
    sdb = 10*log10(Pk_lin(:,k) + 1e-30);
    sdb = sdb - max(sdb);
    plot(fMHz, sdb, 'Color', colors_odd(i,:), 'LineWidth', 1.8);
end
xline(-BW/2/1e6, 'k:', 'LineWidth', 1);
xline( BW/2/1e6, 'k:', 'LineWidth', 1);
xline(-3*BW/2/1e6, ':', 'Color',[.55 .55 .55]);
xline( 3*BW/2/1e6, ':', 'Color',[.55 .55 .55]);
xline(-5*BW/2/1e6, ':', 'Color',[.75 .75 .75]);
xline( 5*BW/2/1e6, ':', 'Color',[.75 .75 .75]);
xlabel('Frequency (MHz)');
ylabel('Normalized PSD (dB, each curve to own peak)');
title('Per-Order Spectrum -- Odd Orders:  y_k = x \cdot |x|^{k-1}');
legend(arrayfun(@(kk)sprintf('k=%d',kk), odd_orders,'UniformOutput',false), ...
       'Location','northeast');
xlim([-120 120]); ylim([-100 5]); grid on;
text(11, 2, 'ACLR1', 'FontSize', 9, 'Color', [.4 .4 .4]);
text(31, 2, 'ACLR2', 'FontSize', 9, 'Color', [.4 .4 .4]);
saveas(gcf, 'PA_poly_odd.png');
fprintf('Saved PA_poly_odd.png\n');

%% 7. Figure 2 -- All-Order Overlay (normalized, odd solid / even dashed)
all_orders = 1:9;
colors_all = [
    0.85 0.10 0.10;    % k=1 red  (odd)
    0.20 0.40 0.90;    % k=2 blue (even)
    0.95 0.45 0.10;    % k=3 orange (odd)
    0.20 0.70 0.80;    % k=4 cyan (even)
    0.85 0.70 0.10;    % k=5 yellow-gold (odd)
    0.30 0.60 0.30;    % k=6 green (even)
    0.65 0.25 0.10;    % k=7 dark orange (odd)
    0.40 0.80 0.60;    % k=8 teal (even)
    0.45 0.15 0.10];   % k=9 brown (odd)
linestyles = {'-','--','-','--','-','--','-','--','-'};

figure('Position',[80 80 1250 650]);
hold on;
for k = all_orders
    sdb = 10*log10(Pk_lin(:,k) + 1e-30);
    sdb = sdb - max(sdb);
    plot(fMHz, sdb, 'Color', colors_all(k,:), 'LineWidth', 1.6, ...
         'LineStyle', linestyles{k});
end
xline(-BW/2/1e6, 'k:', 'LineWidth',1);
xline( BW/2/1e6, 'k:', 'LineWidth',1);
xline(-3*BW/2/1e6, ':', 'Color',[.55 .55 .55]);
xline( 3*BW/2/1e6, ':', 'Color',[.55 .55 .55]);
xline(-5*BW/2/1e6, ':', 'Color',[.75 .75 .75]);
xline( 5*BW/2/1e6, ':', 'Color',[.75 .75 .75]);
xlabel('Frequency (MHz)');
ylabel('Normalized PSD (dB, each curve to own peak)');
title('Per-Order Spectrum -- All Orders  (odd = solid, even = dashed)');
lbl = cell(1, numel(all_orders));
for k = all_orders
    if mod(k,2)==0
        lbl{k} = sprintf('k=%d (even)', k);
    else
        lbl{k} = sprintf('k=%d', k);
    end
end
legend(lbl, 'Location','eastoutside');
xlim([-120 120]); ylim([-100 5]); grid on;
saveas(gcf, 'PA_poly_all.png');
fprintf('Saved PA_poly_all.png\n');

%% 8. Figure 3 -- Per-Order 3x3 Grid (absolute dB, reference = k=1 peak)
ref_peak = max(Pk_lin(:,1));

figure('Position',[50 50 1400 900]);
for k = all_orders
    subplot(3,3,k);
    yl = [-120 10];

    % shaded in-band region
    patch([-BW/2/1e6 BW/2/1e6 BW/2/1e6 -BW/2/1e6], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.90 0.95 1.0], ...
          'EdgeColor','none'); hold on;

    sdb = 10*log10(Pk_lin(:,k)/ref_peak + 1e-30);
    plot(fMHz, sdb, 'Color', [0.10 0.30 0.80], 'LineWidth', 1.3);

    xline(-BW/2/1e6, 'k--', 'LineWidth', 0.8);
    xline( BW/2/1e6, 'k--', 'LineWidth', 0.8);
    xline(-3*BW/2/1e6, ':', 'Color',[.55 .55 .55]);
    xline( 3*BW/2/1e6, ':', 'Color',[.55 .55 .55]);

    if mod(k,2) == 0
        par_str = ' (even)';
    else
        par_str = '';
    end
    title(sprintf('k = %d%s   (expected BW \\approx %d MHz)', ...
                  k, par_str, k*BW/1e6), 'FontSize', 10);
    xlabel('Frequency (MHz)', 'FontSize', 9);
    ylabel('PSD (dB, ref=k=1 peak)', 'FontSize', 9);
    xlim([-120 120]); ylim(yl); grid on;
end
sgtitle('Per-Order Spectral Regrowth  (y_k = x \cdot |x|^{k-1},  a_k = 1)');
saveas(gcf, 'PA_poly_grid.png');
fprintf('Saved PA_poly_grid.png\n');

%% 9. ACLR Table
inband_mask = abs(f) <= BW/2;
acl1_mask   = (abs(f) >  BW/2)   & (abs(f) <= 3*BW/2);
acl2_mask   = (abs(f) >  3*BW/2) & (abs(f) <= 5*BW/2);

P1_in = sum(Pk_lin(inband_mask, 1));      % reference: k=1 in-band power

fprintf('\n=== Per-Order Power Table  (a_k = 1) ===\n');
fprintf('In-band = |f| <= 10 MHz ;  ACLR1 = 10..30 MHz ;  ACLR2 = 30..50 MHz\n');
fprintf('P_in(dB) is referenced to k=1 in-band power.\n\n');
fprintf('%5s %12s %12s %12s %15s\n', ...
        'Order', 'P_in (dB)', 'ACLR1(dBc)', 'ACLR2(dBc)', 'Expected BW');
fprintf('%5s %12s %12s %12s %15s\n', ...
        '-----', '---------', '---------', '---------', '-----------');
for k = all_orders
    Pin = sum(Pk_lin(inband_mask, k));
    Pa1 = sum(Pk_lin(acl1_mask,   k));
    Pa2 = sum(Pk_lin(acl2_mask,   k));

    Pin_db = 10*log10(Pin / P1_in);
    if Pa1 < 1e-25
        aclr1 = -Inf;
    else
        aclr1 = 10*log10(Pa1 / Pin);
    end
    if Pa2 < 1e-25
        aclr2 = -Inf;
    else
        aclr2 = 10*log10(Pa2 / Pin);
    end

    fprintf('%5d %12.2f %12.2f %12.2f %12d MHz\n', ...
            k, Pin_db, aclr1, aclr2, k*BW/1e6);
end

%% 10. Even-Order Contribution to ACLR -- Discussion
fprintf('\n=== Even-Order Contribution to ACLR ===\n\n');

fprintf('[Physical RF (bandpass) interpretation]\n');
fprintf('  At RF, y(t) = sum a_k x(t)^k is a real polynomial of a real\n');
fprintf('  passband signal.  For ODD k, the IM products land near the\n');
fprintf('  carrier fc -> observable in-band and ACLR distortion.  For\n');
fprintf('  EVEN k, the products land at DC and 2*fc, which are filtered\n');
fprintf('  out by the PA output matching network / antenna bandpass.\n');
fprintf('  At the antenna, only ODD orders produce measurable ACLR.\n');
fprintf('  --> User intuition CORRECT.\n\n');

fprintf('[Baseband polynomial (this simulation) interpretation]\n');
fprintf('  In the complex baseband model y_k = x * |x|^(k-1), even-k\n');
fprintf('  terms contain |x|^(odd), a non-polynomial (sqrt-based)\n');
fprintf('  envelope factor that DOES produce energy centered at DC.\n');
fprintf('  The ACLR1 / ACLR2 columns for k = 2, 4, 6, 8 in the table\n');
fprintf('  above are NOT -Inf.  This is a mathematical artifact of the\n');
fprintf('  baseband polynomial basis, not a physical effect.\n\n');

fprintf('[Practical implication for DPD]\n');
fprintf('  DPD polynomial basis sets conventionally use ONLY odd orders\n');
fprintf('  (k = 1, 3, 5, ...).  Reason: the physical PA already filters\n');
fprintf('  out even-order products, so including them in the DPD model\n');
fprintf('  is over-parameterization -- the even coefficients end up\n');
fprintf('  fitting noise/artifacts, not actual PA behavior.\n\n');

fprintf('[Where each order''s regrowth lands]\n');
fprintf('  Order k basis y_k = x*|x|^(k-1) has baseband bandwidth ~ k*BW.\n');
fprintf('  So energy boundaries (one-sided):\n');
for k = all_orders
    fprintf('    k=%d  : energy within  |f| < %3d MHz  (%.0f MHz two-sided)\n', ...
            k, k*BW/2/1e6, k*BW/1e6);
end
fprintf('\n  -> ACLR1 (10..30 MHz) is populated starting from k=2 / k=3.\n');
fprintf('  -> ACLR2 (30..50 MHz) is populated starting from k=4 / k=5.\n');
fprintf('  -> Orders k >= 7 have most of their energy outside ACLR2.\n');

fprintf('\nAll figures saved to current directory.  Done.\n');
