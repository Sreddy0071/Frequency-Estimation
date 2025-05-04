Eigen SNR variance
N = 100; % Define the time vector
n = 0:N-1;
% Signal generation: sum of complex exponentials in noise
frequencies = [0.25, 0.5, 0.75]; % Normalized frequencies (x pi rad/sample)
amplitudes = [1.5, 2, 1.5];
phases = [0, pi/4, pi/2];
signal = amplitudes(1)*exp(1i*2*pi*frequencies(1)*n + phases(1)) + amplitudes(2)*exp(1i*2*pi*frequencies(2)*n + phases(2)) + amplitudes(3)*exp(1i*2*pi*frequencies(3)*n + phases(3));
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
%MSE vs SNR
% Parameters
numSignals = 2; % Number of signals
numAntennas = 8; % Number of antennas in the array
numSamples = 1024; % Number of samples
snrValues = -20:5:20; % SNR range in dB
trueDOAs = [30, 60]; % True directions of arrival in degrees
wavelength = 1; % Wavelength and array element spacing
d = wavelength / 2;
angles = -90:1:90; % Angles for spectrum
mseMUSIC = zeros(length(snrValues), 1); % Initialize MSE storage
mseEIG = zeros(length(snrValues), 1);
for idx = 1:length(snrValues) % Main simulation
fprintf('Processing SNR = %d dB\n', snrValues(idx)); % Diagnostic print
% Generate signals and noise
noisePower = 10^(-snrValues(idx) / 10);
noise = sqrt(noisePower / 2) * (randn(numAntennas, numSamples) + 1i * randn(numAntennas, numSamples));
signals = zeros(numAntennas, numSamples);
for k = 1:numSignals
angleRad = trueDOAs(k) * pi / 180;
steeringVector = exp(1i * 2 * pi * d / wavelength * sin(angleRad) * (0:numAntennas-1).');
signals = signals + steeringVector * (randn(1, numSamples) + 1i * randn(1, numSamples));
end
receivedSignal = signals + noise;
R = (receivedSignal * receivedSignal') / numSamples; % Spatial correlation matrix
% Eigendecomposition
[eigenvectors, eigenvalues] = eig(R);
[eigenvalues, ind] = sort(diag(eigenvalues), 'descend');
eigenvectors = eigenvectors(:, ind);
signalSubspace = eigenvectors(:, 1:numSignals);
noiseSubspace = eigenvectors(:, numSignals+1:end);
% MUSIC Spectrum
musicSpectrum = zeros(size(angles));
for ang = 1:length(angles)
angleRad = angles(ang) * pi / 180;
steeringVector = exp(1i * 2 * pi * d / wavelength * sin(angleRad) * (0:numAntennas-1).');
musicSpectrum(ang) = 1 / (steeringVector' * noiseSubspace * noiseSubspace' * steeringVector);
end
musicSpectrum = 10 * log10(abs(musicSpectrum) / max(abs(musicSpectrum)));
% Eigenvector Method Spectrum
eigSpectrum = zeros(size(angles));
for ang = 1:length(angles)
angleRad = angles(ang) * pi / 180;
steeringVector = exp(1i * 2 * pi * d / wavelength * sin(angleRad) * (0:numAntennas-1).');
eigSpectrum(ang) = abs(steeringVector' * signalSubspace * signalSubspace' * steeringVector);
end
eigSpectrum = 10 * log10(eigSpectrum / max(eigSpectrum));
% Find peaks and compute MSE
[~, doaIndicesMUSIC] = findpeaks(musicSpectrum, 'SortStr', 'descend', 'NPeaks',numSignals);
estimatedDOAsMUSIC = angles(doaIndicesMUSIC);
mseMUSIC(idx) = mean((sort(estimatedDOAsMUSIC) - sort(trueDOAs)).^2);
[~, doaIndicesEIG] = findpeaks(eigSpectrum, 'SortStr', 'descend', 'NPeaks',numSignals);
estimatedDOAsEIG = angles(doaIndicesEIG);
mseEIG(idx) = mean((sort(estimatedDOAsEIG) - sort(trueDOAs)).^2);
end
% Plotting MSE vs SNR for MUSIC
figure;
plot(snrValues, 10 * log10(mseMUSIC + 1e-10), 'b-o', 'LineWidth', 2, 'MarkerSize', 8,'DisplayName', 'MUSIC Method');
title('MSE vs. SNR for MUSIC Method');
xlabel('SNR (dB)');
ylabel('MSE (dB)');
legend show;
% Plotting MSE vs SNR for Eigenvector Method
figure;
plot(snrValues, 10 * log10(mseEIG + 1e-10), 'r-s', 'LineWidth', 2, 'MarkerSize', 8,'DisplayName', 'Eigenvector Method');
title('MSE vs. SNR for Eigenvector Method');
xlabel('SNR (dB)');
ylabel('MSE (dB)');
legend show;