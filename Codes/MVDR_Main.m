
    clear;
    clc;
    close all;

    % Parameters
    N = 10;                         % number of antennas
    noise = 3;                      % Noise level
    P = 10;                         % number of sources
    lambda = 5;                     % wavelength of the signal
    d = lambda/2;                   % distance between the antennas

    % SNR range
    snr_range = -20:5:30;

    % SNR values
    snr_values = -10:5:30;

    % Signal Parameters
    sig = [10^2, 10^1.5, 10^2.2, 10^1.8, 10^1.7, 10^1.9, 10^1.6, 10^1.4, 10^1.3, 10^1.2];
    theta = [-10, 20, 45, -60, 70, 0, 50, -40, -80, 90];

    % Array position
    position = zeros(2, N);
    position(2,:) = 0;
     for mm = 1:N
    position(1,mm) = (mm-1)*lambda./2;      % Element positions in the array
     end

  % covariance matrix with noise
   Rx0 = noise*eye(N);

  % Signals
for i = 1:10
    v = [-sind(theta(i)), cosd(theta(i))];
    k = 2*pi/lambda*v;
    as = exp(-1i*(k*position));
    as = as.';
    Rx0 = Rx0 + sig(i)*as*as';
end

    % Define theta_0
    theta_0 = -90:1:90;

    % MVDR beamforming response
    P_mvdr = mvdr_beamforming_response(Rx0, position, lambda, theta_0);


    % RMSE for different SNR values
    rmse_mvdr = zeros(size(snr_values));
    rmse_symbols = zeros(size(snr_values));              % RMSE of recovered symbols
    

   for i = 1:length(snr_values)
    snr = snr_values(i);
    
    % Adding noise to the received signal covariance matrix
    noise_var = 10^(-snr/10);
    Rx0_noisy = Rx0 + noise_var*eye(N);
    
    % MVDR beamforming response
    P_mvdr_snr = mvdr_beamforming_response(Rx0_noisy, position, lambda, theta_0);
    rmse_mvdr(i) = sqrt(mean((abs(P_mvdr_snr) - abs(P_mvdr)).^2));
    
    
    % Plot MVDR beamforming response
    plot(theta_0, 10*log10(abs(P_mvdr_snr)), 'LineWidth', 1.5); hold on;
 
   % QPSK Modulation
    symbols = randi([0, 3], 1, N);                                  %  symbols for SNR
    modulatedSymbols = exp(1i*2*pi*symbols/4);

   % Transmitter 
   transmittedSymbols = modulatedSymbols .* mean(abs(P_mvdr_snr));     % Modulating with MVDR

   % Adding Gaussian noise
    receivedSymbols = transmittedSymbols + sqrt(noise_var/2)*(randn(size(transmittedSymbols))+1i*randn(size(transmittedSymbols)));

    % Receiver (Demodulation)
    demodulatedSymbols = round(angle(receivedSymbols)/(2*pi)*4);
    demodulatedSymbols(demodulatedSymbols<0) = demodulatedSymbols(demodulatedSymbols<0) + 4;

    % BER Calculation
    numBitErrors = sum(abs(demodulatedSymbols - symbols));
    ber = numBitErrors / N;
    disp(['Bit error rate = ' num2str(ber)]);

    % Calculate RMSE of recovered symbols
    rmse_symbols(i) = sqrt(mean((demodulatedSymbols - symbols).^2));
    
    % Check symbols match
    numErrors = sum(symbols ~= demodulatedSymbols);
    disp(['Number of errors = ', num2str(numErrors)]);

   end

   % Plot
   title('MVDR Beamforming Response at Different SNRs');
   xlabel('Theta (degrees)');
   ylabel('Response (dB)');
   legend(arrayfun(@(snr) sprintf('SNR = %d dB', snr), snr_values, 'UniformOutput', false));
   grid on;
   
   
    % Print RMSE values
    disp("RMSE values (MVDR):");
    disp(rmse_mvdr);

    % Plot RMSE of DOA vs SNR
    figure;
    plot(snr_values, rmse_mvdr, 'LineWidth', 1.5);
    title('RMSE of DOA versus SNR');
    xlabel('SNR (dB)');
    ylabel('RMSE (dB)');
    grid on;

    % Plot RMSE of recovered symbols vs SNR
    figure;
    plot(snr_values, rmse_symbols, 'LineWidth', 1.5, 'LineStyle', '--');
    title('RMSE of Recoverd Symbols versus SNR');
    xlabel('SNR (dB)');
    ylabel('RMSE (dB)');
    grid on;

  % RMSE arrays for each source and peak arrays
   peak_capon_array = zeros(1, length(snr_range));

for snr_idx = 1:length(snr_range)
    
   %  Array steering matrix
    A = exp(-1i*2*pi*d*(0:N-1)'*sind(theta)/lambda);

    % Generate the noise signal 
    snr = snr_range(snr_idx);
    noise = ones(N,1000)*sqrt(10^(-snr/10));

    % Generate the signals received 
    received_signal_1 = A * exp(1i*2*pi*rand(P, 1000));
  
    % Add the noise to the signals
    received_signal_1 = received_signal_1 + noise;
 
    % Beamforming
    theta_scan = -180:0.1:180;                                           % scanning angles
    A_scan = exp(-1i*2*pi*d*(0:N-1)'*sind(theta_scan)/lambda);           % steering matrix for scanning angles

   % Capon Beamforming and peak search
   capon_output = capon_beamforming_method(received_signal_1, A_scan);  
   [peak_capon, angle_capon] = max(capon_output(:));
   peak_capon_array(snr_idx) = peak_capon;   
   peak_capon_array = round(abs(peak_capon_array));                       % Store peak for Capon method                           
end

    % Print peak value
   disp("Peak value (Capon):");
   disp(peak_capon_array);

% Plot Capon beamforming response and peak
figure;
plot(theta_scan, 10*log10(abs(capon_output)), 'LineWidth', 1.5);
hold on;
plot(theta_scan(angle_capon), 10*log10(abs(peak_capon)), 'r*', 'MarkerSize', 10);
title('Capon Beamforming Response and Peak');
xlabel('Theta (degrees)');
ylabel('Response (dB)');
legend('Capon Response', 'Peak');
grid on;

%Plot
figure;
scatter(real(receivedSymbols), imag(receivedSymbols));    % Scatter plot of the received symbols
title('Received Symbols in Complex Plane');
xlabel('Real Part');
ylabel('Imaginary Part');
grid on;


% MVDR beamforming response function
function P = mvdr_beamforming_response(Rx0, position, lambda, theta_0)
    N = size(position, 2);
    P = zeros(1,length(theta_0));
    a = zeros(N,1);

    for mm = 1:length(theta_0)
        vd = [-sind(theta_0(mm)), cosd(theta_0(mm))];       % Steering vector angle
        kd = 2*pi./lambda*vd;                               % Wavenumber for the angle
        for nn = 1:N
            a(nn) = exp(-1i*(kd*position(:,nn)));           % Steering vector for each element
        end
        w_beamforming = inv(Rx0)*a/(a'*inv(Rx0)*a);         % Weight 

        P(mm) = w_beamforming'*Rx0*w_beamforming;           % Beamforming response
    end
end



% Implement Capon beamforming method
function capon_output = capon_beamforming_method(received_signal, A_scan)
    
    R = received_signal * received_signal' / size(received_signal,2);  % Covariance matrix
    invR = inv(R + 1e-8*eye(size(R)));                                 
    capon_output = zeros(size(A_scan,2),1);
    for i = 1:size(A_scan,2)
        a = A_scan(:,i);
        capon_output(i) = 1 / (a' * invR * a);
    end
end