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
lambda = 1;                 % wavelength
d = lambda/2;               % distance between the antennas

% SNR range
snr_range = -20:1:20; 

% RMSE arrays for each source
rmse_1 = zeros(1, length(snr_range));
rmse_2 = zeros(1, length(snr_range));

% RMSE array for the symbols
rmse_symbols = zeros(1, length(snr_range));

%BER array
ber = zeros(1, length(snr_range));

% Number of iterations
num_iterations = 50;

% Simulations
num_simulations = 100;

% symbols 
num_symbols = num_iterations;
symbols = randi([0 1], num_symbols, 1);

% BPSK Modulation
modulated_symbols = 2*symbols-1;

for snr_idx = 1:length(snr_range)
    rmse_1_temp = zeros(1, num_simulations);
    rmse_2_temp = zeros(1, num_simulations);
    ber_temp = zeros(1, num_simulations);

    for sim_idx = 1:num_simulations
        
        %  Array steering matrix
        A = exp(-1i*2*pi*d*(0:N-1)'*sind(theta)/lambda);

        % Generating the noise signal
        snr = snr_range(snr_idx);
        noise = ones(N,num_iterations)*sqrt(10^(-snr/10));

        % Recived source
        received_signal_1 = A * (exp(1i*2*pi*rand(P, num_iterations)).*modulated_symbols');
        received_signal_2 = A * (exp(1i*2*pi*rand(P, num_iterations)).*modulated_symbols');


        % Adding noise 
        received_signal_1 = received_signal_1 + noise;
        received_signal_2 = received_signal_2 + noise;

        % Covariance matrix of the received signals
        R_1 = received_signal_1 * received_signal_1' / size(received_signal_1,2);
        R_2 = received_signal_2 * received_signal_2' / size(received_signal_2,2);

        % Use ESPRIT method 
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

        % Calculating RMSE 
        rmse_1_temp(sim_idx) = sqrt(mean((sort(theta_hat_1) - theta_1).^2));
        rmse_2_temp(sim_idx) = sqrt(mean((sort(theta_hat_2) - theta_2).^2));

        % Adding noise
        noisy_symbols = awgn(modulated_symbols, snr, 'measured');
    
        % BPSK Demodulation
        demodulated_symbols = real(noisy_symbols) > 0;
    
        % RMSE for the symbols
        rmse_symbols(snr_idx) = sqrt(mean((symbols - demodulated_symbols).^2));

        % Compute Bit Error Rate (BER)
        ber_temp(sim_idx) = sum(abs(symbols - demodulated_symbols)) / length(symbols);
    end

    % Average RMSE 
    rmse_1(snr_idx) = mean(rmse_1_temp);
    rmse_2(snr_idx) = mean(rmse_2_temp);

    % Average BER 
    ber(snr_idx) = mean(ber_temp);
    
    % Check symbols match
   numErrors = sum(symbols ~= demodulated_symbols);
   fprintf('Number of errors = %d\n', numErrors);
end

% Calculate the combined RMSE
rmse = rmse_1 + rmse_2;

% Print the RMSE values
disp('RMSE for Signal 1:');
disp(rmse_1);
disp('RMSE for Signal 2:');
disp(rmse_2);
disp('Combined RMSE:');
disp(rmse);
disp('RMSE for Symbols:');
disp(rmse_symbols);
disp('Bit Error Rate (BER):');
disp(ber);

% Beamforming using Delay-and-Sum
weights = ones(N, 1);                                                                 
delay_sum_output_1 = abs(sum(bsxfun(@times, received_signal_1, weights), 1));         
delay_sum_output_2 = abs(sum(bsxfun(@times, received_signal_2, weights), 1)); 

% Combine both signals
combined_signal = delay_sum_output_1 + delay_sum_output_2;

%  Peak search 
peak_values = [];
peak_indices = [];

for i = 2:length(combined_signal)-1                 %Single dataset
    if combined_signal(i) > combined_signal(i-1) && combined_signal(i) > combined_signal(i+1)
        peak_values = [peak_values, combined_signal(i)]; 
        peak_indices = [peak_indices, i];
    end
end

% Highest peak and index
[max_peak_value, max_peak_index] = max(peak_values);
 max_peak_location = peak_indices(max_peak_index);

% Display Highest Peak value
disp('Highest Peak Value:');
disp(max_peak_value);
disp('Highest Peak Location:');
disp(max_peak_location);

% Plot RMSE of DOA versus SNR
figure;
plot(snr_range, rmse_1, '-', snr_range, rmse_2, '-', snr_range, rmse, '-');
grid on;
title('RMSE of DOA vs. SNR');
xlabel('SNR (dB)');
ylabel('RMSE');
legend('RMSE (Signal 1)', 'RMSE (Signal 2)', 'Combined RMSE');

% Plot the delay and sum beamforming outputs
figure;
plot(1:size(received_signal_1, 2), delay_sum_output_1);
hold on;
plot(1:size(received_signal_2, 2), delay_sum_output_2);
xlabel('Snapshot');
ylabel('Magnitude');
title('Delay and Sum Beamforming Output for Signal 1 and Signal 2');
legend('Signal 1', 'Signal 2');
grid on;

% Display BER value
disp('Bit Error Rate (BER) for each SNR level:');
disp(ber);


% Plot RMSE of Symbols versus SNR
figure;
plot(snr_range, rmse_symbols, '-');
grid on;
title('RMSE of Recoverd Symbols vs. SNR');
xlabel('SNR (dB)');
ylabel('RMSE');
legend('RMSE Symbols');


% Plot the combined delay and sum beamforming output with highest peak
figure;
plot(1:size(combined_signal, 2), combined_signal);
hold on;
plot(max_peak_location, max_peak_value, 'r*', 'MarkerSize', 10);
xlabel('Snapshot');
ylabel('Magnitude');
title('Combined Delay and Sum Beamforming Output with Highest Peak');
legend('Combined Signal', 'Highest Peak');
grid on;
