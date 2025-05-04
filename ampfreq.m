Music Amp freq
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
% Prepare a figure for the combined plot
figure;
hold on; % Enable plotting on the same axes
for i = 1:size(frequency_sets, 1)
% Generate the signal with the current set of frequencies
frequencies = frequency_sets(i, :);
signal = amplitudes(1)*exp(1i*2*pi*frequencies(1)*n + phases(1)) + amplitudes(2)*exp(1i*2*pi*frequencies(2)*n + phases(2)) + amplitudes(3)*exp(1i*2*pi*frequencies(3)*n + phases(3));
signal = awgn(signal, SNR, 'measured'); % Add white Gaussian noise
% Correlation matrix estimation using the standard covariance method
R = corrmtx(signal, order, 'cov');
% MUSIC algorithm
[Pxx, F] = pmusic(R, order, [], 1); % Spectrum estimation
% Plotting on the same axes
plot(F, 10*log10(Pxx), 'LineWidth', 1.5);
end
hold off; % Disable plotting on the same axes
title('Combined Music Algorithm Spectral Density for Different Frequency Sets');
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Spectral Power Density (dB)');
legend('Frequencies [0.25, 0.5, 0.75] \times \pi', ...
'Frequencies [0.3, 0.6, 0.9] \times \pi', ...
'Frequencies [0.2, 0.4, 0.6] \times \pi');
grid on;