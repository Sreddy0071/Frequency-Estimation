Eigen
% Time vector
n = 0:99;
% Signal generation: Complex sinusoids with noise
s = exp(1i*pi/2*n) + 2*exp(1i*pi/4*n) + exp(1i*pi/3*n) + randn(1,100);
% Generate the correlation matrix using the modified covariance method
X = corrmtx(s, 12, 'mod');
% Spectral estimation using MUSIC algorithm
[pxx, f] = peig(X, 3, 'whole', 512); % Ensure 512 FFT points for higher resolution
% Plotting the spectral density
figure;
plot(f, 10*log10(pxx));
grid on;
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Power Spectral Density (dB)');
title('Power Spectral Density Estimate Using Eigen vector Method');

%Eigen_Amp_variance
% Define the time vector
N = 100;
n = 0:N-1;
% Signal generation parameters
frequencies = [0.25, 0.5, 0.75]; % Normalized frequencies (x pi rad/sample)
phases = [0, pi/4, pi/2];
SNR = 10; % Signal-to-noise ratio in dB
% Varying amplitudes
amplitudeSets = [1, 2, 3; 2, 3, 4; 3, 4, 5];
% Prepare a figure for the combined plot
figure;
hold on; % Enable plotting on the same axes
for i = 1:size(amplitudeSets, 1)
amplitudes = amplitudeSets(i, :);
signal = amplitudes(1)*exp(1i*2*pi*frequencies(1)*n + phases(1)) + amplitudes(2)*exp(1i*2*pi*frequencies(2)*n + phases(2)) + amplitudes(3)*exp(1i*2*pi*frequencies(3)*n + phases(3));
signal = awgn(signal, SNR, 'measured'); % Add white Gaussian noise
% Correlation matrix estimation and Eigenvector algorithm
R = corrmtx(signal, 3, 'cov');
[Pxx, F] = peig(R, 3, 512, 1);
% Plotting on the same axes
plot(F, 10*log10(Pxx), 'LineWidth', 1.5);
end
hold off; % Disable plotting on the same axes

%Eigenfreq_variance
% Define the time vector
N = 100;
n = 0:N-1;
% Initial parameters
amplitudes = [1.5, 2, 1.5];
phases = [0, pi/4, pi/2];
SNR = 10; % Signal-to-noise ratio in dB
order = 3; % Model order, generally set to the number of sinusoids
% Different sets of frequencies to be tested
frequency_sets = [0.25, 0.5, 0.75; 0.3, 0.6, 0.9; 0.2, 0.4, 0.6];
% Set up the figure
figure;
hold on; % Enable plotting on the same axes
for i = 1:size(frequency_sets, 1)
% Generate the signal with current set of frequencies
frequencies = frequency_sets(i, :);
signal = amplitudes(1)*exp(1i*2*pi*frequencies(1)*n + phases(1)) + ...
amplitudes(2)*exp(1i*2*pi*frequencies(2)*n + phases(2)) + ...
amplitudes(3)*exp(1i*2*pi*frequencies(3)*n + phases(3));
signal = awgn(signal, SNR, 'measured'); % Add white Gaussian noise
% Correlation matrix estimation using the standard covariance method
R = corrmtx(signal, order, 'cov');
% MUSIC algorithm
[Pxx, F] = peig(R, order, [], 1); % Spectrum estimation
% Plotting the spectral density
plot(F, 10*log10(Pxx), 'LineWidth', 1.5);
end
hold off; % Disable plotting on the same axes
title('Combined Spectral Density Plots');
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Spectral Power Density (dB)');
legend('Frequencies [0.25, 0.5, 0.75] \times \pi', ...
'Frequencies [0.3, 0.6, 0.9] \times \pi', ...
'Frequencies [0.2, 0.4, 0.6] \times \pi');
grid on;
%Eigen_MSMESNR
% Define the time vector
N = 100;
n = 0:N-1;
% Signal generation: sum of complex exponentials in noise
frequencies = [0.25, 0.5, 0.75]; % Normalized frequencies (x pi rad/sample)
amplitudes = [1.5, 2, 1.5];
phases = [0, pi/4, pi/2];
signal = amplitudes(1)*exp(1i*2*pi*frequencies(1)*n + phases(1)) + ...
amplitudes(2)*exp(1i*2*pi*frequencies(2)*n + phases(2)) + ...
amplitudes(3)*exp(1i*2*pi*frequencies(3)*n + phases(3));
SNR = 10; % Signal-to-noise ratio in dB
signal = awgn(signal, SNR, 'measured'); % Add white Gaussian noise
% Correlation matrix estimation using the standard covariance method
order = 3; % Model order, generally set to the number of sinusoids
R = corrmtx(signal, order, 'cov');
% MUSIC algorithm
[Pxx, F] = peig(R, order, [], 1); % Spectrum estimation
% Plotting the spectral density
figure;
plot(F, 10*log10(Pxx), 'b', 'LineWidth', 1.5);
title('Spectral Density Estimate Using Eigen vector Algorithm');
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Spectral Power Density (dB)');
grid on;