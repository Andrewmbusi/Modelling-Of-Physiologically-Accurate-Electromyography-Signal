% MODELLING OF PHYSIOLOGICALLY ACCURATE EMG SIGNAL
clc;clear;close all;
%% Parameters and Setup
Fs = 2000;              
T = 40;                 
totalSamples = Fs * T;  
time = (0:totalSamples-1) / Fs;
% Bandpass filter (10–400 Hz)
[b, a] = butter(4, [10 400] / (Fs/2), 'bandpass');

%% Generate Physiological EMG
physio_EMG_raw = zeros(1, totalSamples + 5000);
num_MUAPs = 4000;
fprintf('Simulating %d Motor Units...\n', num_MUAPs);
for k = 1:num_MUAPs
    dist = 0.5 + (1.5) * rand();
    muap = generate_MUAP_simple(dist);
    start_idx = randi([1, totalSamples]);
    end_idx = start_idx + length(muap) - 1;
    physio_EMG_raw(start_idx:end_idx) = ...
        physio_EMG_raw(start_idx:end_idx) + muap;

end
physio_EMG = filter(b, a, physio_EMG_raw(1:totalSamples));
physio_EMG = 2 * physio_EMG / max(abs(physio_EMG));
%% Generate Gaussian EMG
gaussian_raw = randn(1, totalSamples);
gaussian_EMG = filter(b, a, gaussian_raw);
gaussian_EMG = gaussian_EMG * ...
    (std(physio_EMG) / std(gaussian_EMG));
%% Plot Results
create_simple_plot(time, physio_EMG, ...
    'Physiologic EMG', 1, Fs);
create_simple_plot(time, gaussian_EMG, ...
    'Gaussian EMG', 2, Fs);
%% Display Signal Characteristics
fprintf('\n--- Signal Characteristics ---\n');
fprintf('%-25s | %-15s | %-15s\n', ...
    'Metric', 'Physiological', 'Gaussian');
fprintf('%-25s | %-15.4f | %-15.4f\n', ...
    'Skewness', ...
    skewness(physio_EMG), ...
    skewness(gaussian_EMG));
fprintf('%-25s | %-15.4f | %-15.4f\n', ...
    'Kurtosis', ...
    kurtosis(physio_EMG), ...
    kurtosis(gaussian_EMG));
fprintf('%-25s | %-15.4f | %-15.4f\n', ...
    'Mean Frequency (Hz)', ...
    meanfreq(physio_EMG, Fs), ...
    meanfreq(gaussian_EMG, Fs));
fprintf('%-25s | %-15.4f | %-15.4f\n', ...
    'Median Frequency (Hz)', ...
    medfreq(physio_EMG, Fs), ...
    medfreq(gaussian_EMG, Fs));

function create_simple_plot(t, signal, typeStr, figNum, Fs)
    figure(figNum);
    set(gcf, 'Position', [150, 150, 800, 850]);
    tiledlayout(3,1);
    % Time Domain
    nexttile;
    plot(t, signal);
    title(typeStr);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
    % Power Spectrum
    nexttile;
    [pxx, w] = periodogram(signal, ...
        rectwin(length(signal)), length(signal));
    plot(w/pi, 10*log10(pxx));
    title([typeStr, ' Power Spectrum']);
    xlabel('Normalized Frequency');
    ylabel('Power (dB)');
    grid on;
    % Histogram
    nexttile;
    histogram(signal, 100);
    title([typeStr, ' Distribution']);
    xlabel('Value');
    ylabel('Samples');
    grid on;
end
function muap_final = generate_MUAP_simple(distance)
    fiber_len = 500;
    x = 1:fiber_len;
    pos_portion = 1 ./ sqrt(distance.^2 + x.^2);
    padding = zeros(1, 1000);
    pos_padded = [padding, pos_portion, padding];
    shift = 150;
    neg_comp = [zeros(1, shift), -pos_padded];
    pos_comp = [pos_padded, zeros(1, shift)];
    muap_final = pos_comp + neg_comp;
end

