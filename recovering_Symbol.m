% Parameters
N = 10;                     % number of antennas
K = 2;                      % number of sources

theta_1 = 15;               % true DOA of the first source
theta_2 = 40;               % true DOA of the second source
theta = [theta_1 theta_2];  % vector of true DOAs
P = length(theta);          % number of sources
lambda = 1;                 % wavelength of the signal
d = lambda/2;               % distance between the antennas

% Generate the array steering matrix
A = exp(-1i*2*pi*d*(0:N-1)'*sind(theta)/lambda);

% Generate the noise signal with a different noise level
noise_level = 10;          % adjust the noise level as desired
noise = ones(N,1000)*noise_level;

% Generate the signals received by the array
received_signal_1 = A * exp(1i*2*pi*rand(P, 1000));
received_signal_2 = A * exp(1i*2*pi*rand(P, 1000));

% Add the noise to the signals
received_signal_1 = received_signal_1 + noise;
received_signal_2 = received_signal_2 + noise;

% SNR values to test
SNR = -10:5:20;

% Initialize RMSE vectors
rmse_esprit = zeros(1, length(SNR));
rmse_music = zeros(1, length(SNR));
rmse_mvdr = zeros(1, length(SNR));

% Loop over SNR values
for snr_index = 1:length(SNR)
    % Add noise to the received signals
    received_signal_1_noisy = awgn(received_signal_1, SNR(snr_index), 'measured');
    
    % Estimate the covariance matrix of the noisy signals
    R_noisy = received_signal_1_noisy * received_signal_1_noisy' / size(received_signal_1_noisy, 2);
    
    % Use ESPRIT method to estimate the DOAs of the sources
    [U, S, V] = svd(R_noisy);
    U1 = U(:, 1:K);
    V1 = V(:, 1:K);
    F_1 = U1' * V1;
    phi_1 = angle(eig(F_1));
    theta_hat_1 = asin(phi_1 / (2 * pi * d / lambda)) * 180 / pi + theta_1;
    
    % Ensure both vectors are column vectors
    theta_hat_1 = theta_hat_1(:);
    theta = theta(:);
    
    % If the vectors have different sizes, pad the smaller one with NaNs
    if length(theta_hat_1) > length(theta)
        theta = [theta; NaN(length(theta_hat_1)-length(theta), 1)];
    elseif length(theta) > length(theta_hat_1)
        theta_hat_1 = [theta_hat_1; NaN(length(theta)-length(theta_hat_1), 1)];
    end

    % Calculate RMSE for ESPRIT approach
    rmse_esprit(snr_index) = sqrt(nanmean((theta_hat_1 - theta).^2));
    
    % Use MUSIC method to estimate the DOAs of the sources
    [~, ~, V] = svd(R_noisy);
    Vn = V(:, K+1:end);
    theta_scan_music = -90:0.1:90;
    A_music = exp(-1i*2*pi*d*(0:N-1)'*sind(theta_scan_music)/lambda);
    music_spectrum = zeros(size(theta_scan_music));
    for ii = 1:length(theta_scan_music)
        a_theta = A_music(:, ii);
        music_spectrum(ii) = 1 / (a_theta' * Vn * Vn' * a_theta);
    end
    
    % Ensure music_spectrum is real-valued
    music_spectrum = abs(music_spectrum);
    music_spectrum = real(music_spectrum);
    
    [~, peaks] = findpeaks(music_spectrum);
    theta_hat_music = theta_scan_music(peaks);
    
    % Calculate RMSE for MUSIC approach
    theta_hat_music = theta_hat_music(:);
    if length(theta_hat_music) > length(theta)
        theta = [theta; NaN(length(theta_hat_music)-length(theta), 1)];
    elseif length(theta) > length(theta_hat_music)
        theta_hat_music = [theta_hat_music; NaN(length(theta)-length(theta_hat_music), 1)];
    end
    rmse_music(snr_index) = sqrt(nanmean((theta_hat_music - theta).^2));
    
    % Calculate the MVDR beamforming response
    theta_scan_mvdr = -180:0.1:180;                                  % scanning angles
    A_scan_mvdr = exp(-1i*2*pi*d*(0:N-1)'*sind(theta_scan_mvdr)/lambda);   % steering matrix for scanning angles
    R_inv = inv(R_noisy);
    P_mvdr = abs(sum(A_scan_mvdr' * R_inv * received_signal_2, 2)).^2;
    
    % Calculate RMSE for MVDR approach
    [~, idx] = max(P_mvdr);
    if length(idx) > 1  % If multiple peaks, take the first one
        idx = idx(1);
    end
    rmse_mvdr(snr_index) = sqrt(nanmean((theta_scan_mvdr(idx) - theta).^2));
end

% Display RMSE values
fprintf('RMSE values:\n');
fprintf('SNR\t\tESPRIT\t\tMUSIC\t\tMVDR\n');
for snr_index = 1:length(SNR)
    fprintf('%d dB\t%.4f\t\t%.4f\t\t%.4f\n', SNR(snr_index), rmse_esprit(snr_index), rmse_music(snr_index), rmse_mvdr(snr_index));
end

% Plot the RMSE versus SNR
figure;
plot(SNR, rmse_esprit, 'b-o', 'LineWidth', 1.5);
hold on;
plot(SNR, rmse_music, 'g-o', 'LineWidth', 1.5);
plot(SNR, rmse_mvdr, 'm-o', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('RMSE');
title('RMSE of DOA Estimation');
legend('ESPRIT', 'MUSIC', 'MVDR');
grid on;
