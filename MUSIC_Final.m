clear;
clc;
close all;

% Parameters
N = 8;                                        % Number of sensors
K = 3;                                        % Number of sources 
T = 20;                                       % Number of snapshots
R = 100;                                      % SNR
theta = [-50, -30, -10, 0, 10, 30, 50, 70, 80, 90] * pi/180;       % Angle of arrival 
d = 0.3;    % Distance between sources

% SNR range
snr_range = -10:5:30;

% Initialize RMSE array
rmse = zeros(size(snr_range));

% Generate the steering matrix
A = zeros(N, K);
for n = 1:N
    for k = 1:K
        A(n, k) = exp(-1i * 2 * pi * (n-1) * d * sin(theta(k)));
    end
end

% Randomly generated source signals
S = randn(K, T);

for snr_idx = 1:length(snr_range)
    % Current SNR value
    R = snr_range(snr_idx);
    
    % Multiplying source signals and steering matrix
    AS = A * S;
    
    % Add white Gaussian noise to AS with the current SNR
    X = awgn(AS, R);
    
    % Beamforming
    weights = ones(N, 1); % Uniform weights for simplicity
    X_beamformed = sum(bsxfun(@times, X, weights), 1);
    
    % Call the MUSIC function
    [music_spectrum, theta_est] = MUSIC(X, K);
    
    % Calculate RMSE
    [~, idx] = sort(abs(theta_est));
    theta_sorted = theta_est(idx);
    rmse(snr_idx) = sqrt(mean((theta_sorted - theta(1:K)).^2)); % Consider only the first K elements of theta
end

% Plotting
theta_scan_deg = -90:0.05:90;
theta_scan = theta_scan_deg * pi/180;
figure(1);
plot(theta_scan_deg, 10*log10(abs(1./music_spectrum)));
hold on;
grid on;
title('MUSIC Spectrum');
xlabel('Angle (degrees)');
ylabel('Magnitude (dB)');
ylim([-50, 50]);

figure(2);
plot(1:T, real(X_beamformed));
hold on;
grid on;
title('Beamformed Signal');
xlabel('Snapshot');
ylabel('Magnitude');

figure(3);
plot(snr_range, rmse, 'o-');
hold on;
grid on;
title('RMSE of DOA Estimates vs. SNR');
xlabel('SNR (dB)');
ylabel('RMSE');

% Peak search on MUSIC Spectrum plot
[~, peaks_music] = findpeaks(abs(1./music_spectrum));
theta_peaks_music = theta_scan_deg(peaks_music);
figure(1);
hold on;
plot(theta_peaks_music, 10*log10(abs(1./music_spectrum(peaks_music))), 'ro', 'MarkerSize', 5, 'DisplayName', 'Peak');
legend('MUSIC Spectrum', 'Peak');
hold off;

% Peak search on Beamformed Signal plot
[~, peaks_beamformed] = findpeaks(real(X_beamformed));
snapshot_peaks_beamformed = 1:T;
snapshot_peaks_beamformed = snapshot_peaks_beamformed(peaks_beamformed);
figure(2);
hold on;
plot(snapshot_peaks_beamformed, real(X_beamformed(peaks_beamformed)), 'ro', 'MarkerSize', 5, 'DisplayName', 'Peak');
legend('Beamformed Signal', 'Peak');
hold off;

% Peak search on RMSE plot
[~, rmse_peak_idx] = min(rmse);
rmse_peak_snr = snr_range(rmse_peak_idx);
figure(3);
hold on;
plot(rmse_peak_snr, rmse(rmse_peak_idx), 'ro', 'MarkerSize', 5);
text(rmse_peak_snr, rmse(rmse_peak_idx), num2str(rmse(rmse_peak_idx)), 'VerticalAlignment', 'bottom');
hold off;

% Print the RMSE values
disp('RMSE values:');
disp(rmse);

% Print the peak value
%disp(['Peak value of RMSE of DOA Estimates: ' num2str(rmse(rmse_peak_idx))]);



% MUSIC Function
function [music_spectrum, theta_est] = MUSIC(X, K)
    [~, T] = size(X);                           % Number of columns
    Rx = (1/T) * X * X';                         % Covariance matrix
    [V, ~] = eig(Rx);                           % Eigen decomposition of Rx
    V2 = V(:, K+1:end);                         % Estimate noise subspace

    theta_scan_deg = -90:0.05:90;
    theta_scan = theta_scan_deg * pi/180;
    music_spectrum = zeros(size(theta_scan));
    for i = 1:length(theta_scan)
        a = exp(-1i * pi * sin(theta_scan(i)) * (0:1:(size(V2, 1)-1))'); % Create steering column vector
        music_spectrum(i) = 1 / (norm(V2' * a, 2).^2); % MUSIC spatial spectrum
    end

    [~, idx] = sort(music_spectrum, 'descend');
    theta_est = theta_scan(idx(1:K));
end
