clear;
clc;
close all;

% Parameters
N = 8;                                                             % Number of sensors
K = 3;                                                             % Number of sources 
T = 10;                                                            % Number of snapshots
theta = [ -50, 10, 30] * pi/180;                                   % Angle of arrival
d = 0.5;                                                           % Distance between sensors

% SNR range
snr_range = -10:10:20;

% RMSE and BER arrays 
rmse = zeros(size(snr_range));
ber = zeros(K, length(snr_range));                                 % BER array 

% Steering matrix
A = zeros(N, length(theta));
for n = 1:N
    for k = 1:length(theta)
        A(n, k) = exp(-1i * 2 * pi * (n-1) * d * sin(theta(k)));
    end
end

% Generating binary source signals
S = randi([0 1], K, T);

% Modulation (BPSK)
S_mod = 2*S - 1;  % maps 0 to -1 and 1 to +1


% Demodulation and recovery
S_all_sources = S_mod;

% RMSE array for recovered symbols
rmse_symbols = zeros(K, length(snr_range));

% MUSIC function
theta_scan_deg = -90:0.05:90;
theta_scan = theta_scan_deg * pi/180;

%plot
figure;
hold on;
grid on;
title('MUSIC Spectrum');
xlabel('Angle (degrees)');
ylabel('Magnitude (dB)');

% Highest peak
highest_peak_val_all_snr = -Inf;
highest_peak_idx_all_snr = -1;
highest_peak_snr = -1;

for snr_idx = 1:length(snr_range)
    % SNR value
    R = snr_range(snr_idx);
    
    % Source signals and steering matrix
    AS = A * S_mod;
    
    % Adding noise to AS 
    X = awgn(AS, R, 'measured');
    
    % MUSIC function
    [music_spectrum, theta_est] = MUSIC(X, K, theta_scan);
    
    % Calculate RMSE
    [~, idx] = sort(music_spectrum, 'descend');
    theta_sorted = theta_scan(idx(1:K));
    rmse(snr_idx) = sqrt(mean((sort(theta_sorted) - theta).^2));
    
    % Plotting MUSIC Spectrum
    plot(theta_scan_deg, 10*log10(abs(1./music_spectrum)));
    
    % Peak finding
    spectrum_dB_music = 10*log10(abs(1./music_spectrum));
    highest_peak_val_music = -Inf;
    highest_peak_idx_music = -1;
    for idx = 2:(length(spectrum_dB_music)-1)              %Double dataset
        if (spectrum_dB_music(idx) > spectrum_dB_music(idx-1)) && (spectrum_dB_music(idx) > spectrum_dB_music(idx+1))
            if spectrum_dB_music(idx) > highest_peak_val_music
                highest_peak_val_music = spectrum_dB_music(idx);
                highest_peak_idx_music = idx;
            end
        end
    end
    if highest_peak_val_music > highest_peak_val_all_snr
        highest_peak_val_all_snr = highest_peak_val_music;
        highest_peak_idx_all_snr = highest_peak_idx_music;
        highest_peak_snr = R;
    end

    for k = 1:K
        % Beamformed Signal
        beamforming_weights = exp(-1i * 2 * pi * (0:N-1)' * d * sin(theta_est(k)));
        X_beamformed = beamforming_weights' * X;

        % Demodulation
        symbols_recovered = sign(real(X_beamformed));

        % RMSE of recovered symbols
        original_symbols = S_all_sources(k,:);  
        rmse_symbols(k, snr_idx) = sqrt(mean((symbols_recovered - original_symbols).^2));
        
        % BER for source
        bit_errors = sum(original_symbols ~= symbols_recovered);
        total_bits = length(original_symbols);
        ber(k, snr_idx) = bit_errors / total_bits;

        % Print BER
        fprintf('SNR: %d dB, Source: %d, BER: %f\n', R, k, ber(k, snr_idx));
        
        % Check symbols match
        numErrors = sum(original_symbols ~= symbols_recovered);
        fprintf('Number of errors = %d\n', numErrors); 
    end
end

%plot
plot(theta_scan_deg(highest_peak_idx_all_snr), highest_peak_val_all_snr, 'r*', 'MarkerSize', 10, 'DisplayName', 'Peak');
disp(['Highest peak at SNR ', num2str(highest_peak_snr), ' dB, angle ', num2str(theta_scan_deg(highest_peak_idx_all_snr)), ' degrees with value ', num2str(highest_peak_val_all_snr), ' dB']);
xlim([-90, 90]);
legend show;

%plot
figure;
semilogy(snr_range, rmse);
grid on;
title('RMSE of DOA Estimates vs. SNR');
xlabel('SNR (dB)');
ylabel('RMSE');


%plot
figure;
for k = 1:K
    semilogy(snr_range, rmse_symbols(k,:), 'DisplayName', ['Source ', num2str(k)]);
    hold on;
end
grid on;
title('RMSE of Recovered Symbols vs. SNR');
xlabel('SNR (dB)');
ylabel('RMSE');
legend show;


% Plotting BER
figure;
grid on;
hold on;
for k = 1:K
    plot(snr_range, ber(k,:), 'DisplayName', ['Source ', num2str(k)]);
end
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Bit Error Rate (BER) vs SNR');
legend show;



% MUSIC Function
function [music_spectrum, theta_est] = MUSIC(X, K, theta_scan)
    [~, T] = size(X);                                  % Number of columns
    Rx = (1/T) * X * X';                               % Covariance matrix
    [V, ~] = eig(Rx);                                  % Eigen decomposition of Rx
    Vn = V(:, 1:end-K);                                % Estimate noise subspace

    music_spectrum = zeros(size(theta_scan));
    for i = 1:length(theta_scan)
        a = exp(-1i * 2 * pi * (0:1:size(Vn, 1)-1)' * sin(theta_scan(i)));     % Steering column vector
        music_spectrum(i) = 1 / (norm(Vn' * a, 2).^2);                          % MUSIC spatial spectrum
    end

    [~, idx] = sort(music_spectrum, 'descend');
    theta_est = theta_scan(idx(1:K));
end