clear;
clc;
close all;

% Parameters
N = 10;                     % number of antennas
K = 2;                      % number of sources

theta_1 = 15;               % true DOA of the first source
theta_2 = 40;               % true DOA of the second source
theta = [theta_1 theta_2];  % vector of true DOAs
P = length(theta);          % number of sources
lambda = 1;                 % wavelength of the signal
d = lambda/2;               % distance between the antennas

% SNR range
snr_range = -20:5:20;

% Initialize RMSE arrays for each source
rmse_1 = zeros(1, length(snr_range));
rmse_2 = zeros(1, length(snr_range));

for snr_idx = 1:length(snr_range)
    % Generate the array steering matrix
    A = exp(-1i*2*pi*d*(0:N-1)'*sind(theta)/lambda);

    % Generate the noise signal with the current SNR
    snr = snr_range(snr_idx);
    noise = ones(N,1000)*sqrt(10^(-snr/10));

    % Generate the signals received by the array
    received_signal_1 = A * exp(1i*2*pi*rand(P, 1000));
    received_signal_2 = A * exp(1i*2*pi*rand(P, 1000));

    % Add the noise to the signals
    received_signal_1 = received_signal_1 + noise;
    received_signal_2 = received_signal_2 + noise;
    
    % Estimate the covariance matrix of the received signals
    R_1 = received_signal_1 * received_signal_1' / size(received_signal_1,2);
    R_2 = received_signal_2 * received_signal_2' / size(received_signal_2,2);

    % Use ESPRIT method to estimate the DOAs of the sources
    [U,S,V] = svd(R_1);
    U1 = U(:,1:K);
    V1 = V(:,1:K);
    F_1 = U1' * V1;
    phi_1 = angle(eig(F_1));
    theta_hat_1 = asin(phi_1/(2*pi*d/lambda))*180/pi + theta_1;

    [U,S,V] = svd(R_2);
    U1 = U(:,1:K);
    V1 = V(:,1:K);
    F_2 = U1' * V1;
    phi_2 = angle(eig(F_2));
    theta_hat_2 = asin(phi_2/(2*pi*d/lambda))*180/pi + theta_2;

    % Calculate RMSE for each estimated DOA
    rmse_1(snr_idx) = sqrt(mean((sort(theta_hat_1) - theta_1).^2));
    rmse_2(snr_idx) = sqrt(mean((sort(theta_hat_2) - theta_2).^2));
end

% Calculate the combined RMSE
rmse = rmse_1 + rmse_2;

% Print the RMSE values
disp('RMSE values:');
disp(rmse);

% Find the peak value of the full beamforming output
peak_value = max(rmse);

% Print the peak value
disp(['Peak value of full beamforming output: ' num2str(peak_value)]);

% Plot RMSE of DOA versus SNR
figure;
plot(snr_range, rmse_1, 'o-', snr_range, rmse_2, 'o-', snr_range, rmse, 'o-');
grid on;
title('RMSE of DOA vs. SNR');
xlabel('SNR (dB)');
ylabel('RMSE');
legend('RMSE (Signal 1)', 'RMSE (Signal 2)', 'Combined RMSE');



% Beamforming using Delay-and-Sum
weights = ones(N, 1); % Uniform weights for delay and sum beamforming
delay_sum_output_1 = abs(sum(bsxfun(@times, received_signal_1, weights), 1)); % Delay and sum beamforming output for signal 1
delay_sum_output_2 = abs(sum(bsxfun(@times, received_signal_2, weights), 1)); % Delay and sum beamforming output for signal 2

% Find peaks in the delay and sum beamforming outputs
[~, ds_peaks_1] = findpeaks(delay_sum_output_1, 'MinPeakDistance', 1);
[~, ds_peaks_2] = findpeaks(delay_sum_output_2, 'MinPeakDistance', 1);

% Print the peaks of Delay-and-Sum beamforming outputs
disp('Peaks in Delay-and-Sum Beamforming Output (Signal 1):');
disp(ds_peaks_1);
disp('Peaks in Delay-and-Sum Beamforming Output (Signal 2):');
disp(ds_peaks_2);

% Plot the delay and sum beamforming outputs
figure;
plot(1:size(received_signal_1, 2), delay_sum_output_1);
hold on;
plot(1:size(received_signal_2, 2), delay_sum_output_2);
xlabel('Snapshot');
ylabel('Magnitude');
title('Delay and Sum Beamforming Output');
legend('Signal 1', 'Signal 2');
grid on;

% Display the true and estimated DOAs
disp('True DOAs:');
disp(theta);
disp('Estimated DOAs (Signal 1):');
disp(sort(theta_hat_1));
disp('Estimated DOAs (Signal 2):');
disp(sort(theta_hat_2));
