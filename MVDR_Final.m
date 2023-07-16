clear;
clc;
close all;

% Parameters
N = 10;                     % number of antennas
K = 2;                      % number of sources
noise = 3;        % Noise level

theta_1 = 15;               % true DOA of the first source
theta_2 = 40;               % true DOA of the second source
theta = [theta_1 theta_2];  % vector of true DOAs
P = length(theta);          % number of sources
lambda = 1;                 % wavelength of the signal
d = lambda/2;               % distance between the antennas

% SNR range
snr_range = -20:5:20;

% SNR values
snr_values = -10:5:20;

% Signal 1
sig1 = 10^2;
theta1 = -10;
v1 = [-sind(theta1), cosd(theta1)];
k1 = 2*pi/lambda*v1;

% Signal 2
sig2 = 10^1.5;
theta2 = 20;
v2 = [-sind(theta2), cosd(theta2)];
k2 = 2*pi/lambda*v2;


% Signal 3
sig3 = 10^2.2;
theta3 = 45;
v3 = [-sind(theta3), cosd(theta3)];
k3 = 2*pi/lambda*v3;

% Signal 4
sig4 = 10^1.8;
theta4 = -60;
v4 = [-sind(theta4), cosd(theta4)];
k4 = 2*pi/lambda*v4;

% Signal 5
sig5 = 10^1.7;
theta5 = 70;
v5 = [-sind(theta5), cosd(theta5)];
k5 = 2*pi/lambda*v5;

% Signal 6
sig6 = 10^1.9;
theta6 = 0;
v6 = [-sind(theta6), cosd(theta6)];
k6 = 2*pi/lambda*v6;

% Signal 7
sig7 = 10^1.6;
theta7 = 50;
v7 = [-sind(theta7), cosd(theta7)];
k7 = 2*pi/lambda*v7;

% Signal 8
sig8 = 10^1.4;
theta8 = -40;
v8 = [-sind(theta8), cosd(theta8)];
k8 = 2*pi/lambda*v8;

% Signal 9
sig9 = 10^1.3;
theta9 = -80;
v9 = [-sind(theta9), cosd(theta9)];
k9 = 2*pi/lambda*v9;

% Signal 10
sig10 = 10^1.2;
theta10 = 90;
v10 = [-sind(theta10), cosd(theta10)];
k10 = 2*pi/lambda*v10;

% Array position
position = zeros(2, N);
position(2,:) = 0;
for mm = 1:N
    position(1,mm) = (mm-1)*lambda./2;  % Element positions in the array
end

as1 = exp(-1i*(k1*position));   % Signal 1
as2 = exp(-1i*(k2*position));   % Signal 2
as3 = exp(-1i*(k3*position));   % Signal 3
as4 = exp(-1i*(k4*position));   % Signal 4
as5 = exp(-1i*(k5*position));   % Signal 5
as6 = exp(-1i*(k6*position));   % Signal 6
as7 = exp(-1i*(k7*position));   % Signal 7
as8 = exp(-1i*(k8*position));   % Signal 8
as9 = exp(-1i*(k9*position));   % Signal 9
as10 = exp(-1i*(k10*position)); % Signal 10

as1 = as1.';
as2 = as2.';
as3 = as3.';
as4 = as4.';
as5 = as5.';
as6 = as6.';
as7 = as7.';
as8 = as8.';
as9 = as9.';
as10 = as10.';

as = as1 + as2 + as3 + as4 + as5 + as6 + as7 + as8 + as9 + as10;

% Received signal covariance matrix
Rx0 = sig1*as1*as1' + sig2*as2*as2' + sig3*as3*as3' + sig4*as4*as4' + sig5*as5*as5' + ...
      sig6*as6*as6' + sig7*as7*as7' + sig8*as8*as8' + sig9*as9*as9' + sig10*as10*as10' + noise*eye(N);

% Define theta_0
theta_0 = -90:1:90;

% MVDR beamforming response
P_mvdr = mvdr_beamforming_response(Rx0, position, lambda, theta_0);

% Plot MVDR beamforming response
figure
plot(theta_0, 10*log10(abs(P_mvdr)), 'LineWidth', 1.5);
hold on;
plot([theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10], ...
     max(10*log10(abs(P_mvdr)))*ones(1,10), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Angle (degrees)');
ylabel('Power (dB)');
title('MVDR Beamforming Response');
legend('Beamforming Response', 'Signal Angles');

% Compute RMSE for different SNR values
rmse_mvdr = zeros(size(snr_values));

for i = 1:length(snr_values)
    snr = snr_values(i);
    
    % Add noise to the received signal covariance matrix
    noise_var = 10^(-snr/10);
    Rx0_noisy = Rx0 + noise_var*eye(N);
    
    % MVDR beamforming response
    P_mvdr_snr = mvdr_beamforming_response(Rx0_noisy, position, lambda, theta_0);
    rmse_mvdr(i) = sqrt(mean((abs(P_mvdr_snr) - abs(P_mvdr)).^2));
end

% Print RMSE values
disp("RMSE values (MVDR):");
disp(rmse_mvdr);

% Plot RMSE of DOA versus SNR
figure;
plot(snr_values, rmse_mvdr, 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('RMSE (dB)');
title('RMSE of DOA versus SNR');
legend('MVDR');

% Initialize RMSE arrays for each source and peak arrays
rmse_capon = zeros(1, length(snr_range));
peak_capon_array = zeros(1, length(snr_range));

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

    % Beamforming
    theta_scan = -180:0.1:180;  % scanning angles
    A_scan = exp(-1i*2*pi*d*(0:N-1)'*sind(theta_scan)/lambda);   % steering matrix for scanning angles

    % Capon Beamforming and peak search
    capon_output = capon_beamforming_method(received_signal_1, A_scan);  % Implement your actual Capon method here
    [peak_capon, angle_capon] = max(capon_output(:));
    peak_capon_array(snr_idx) = peak_capon;  % Store peak for Capon method
end

% Print peak value
disp("Peak value (Capon):");
disp(peak_capon_array);

% Plotting Peaks of Beamforming Methods vs. SNR
figure;
hold on;
plot(snr_range, real(peak_capon_array), 'o-', 'DisplayName', 'Capon');
grid on;
title('Peaks of Beamforming Methods vs. SNR');
xlabel('SNR (dB)');
ylabel('Peak');

function capon_output = capon_beamforming_method(received_signal, A_scan)
    % Implement Capon beamforming method
    R = received_signal * received_signal' / size(received_signal,2);  % Covariance matrix
    invR = inv(R + 1e-8*eye(size(R)));  % Regularization with small positive value
    capon_output = zeros(size(A_scan,2),1);
    for i = 1:size(A_scan,2)
        a = A_scan(:,i);
        capon_output(i) = 1 / (a' * invR * a);
    end
end

% MVDR beamforming response function
function P = mvdr_beamforming_response(Rx0, position, lambda, theta_0)
    N = size(position, 2);
    P = zeros(1,length(theta_0));
    a = zeros(N,1);

    for mm = 1:length(theta_0)
        vd = [-sind(theta_0(mm)), cosd(theta_0(mm))];  % Steering vector for the desired angle
        kd = 2*pi./lambda*vd;  % Wavenumber for the desired angle
        for nn = 1:N
            a(nn) = exp(-1i*(kd*position(:,nn)));  % Steering vector for each element
        end
        w_beamforming = inv(Rx0)*a/(a'*inv(Rx0)*a);  % Weight 

        P(mm) = w_beamforming'*Rx0*w_beamforming;  % Beamforming response
    end
end
