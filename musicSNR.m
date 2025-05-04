%Music SNR variance
% Signal parameters
N = 100;
n = 0:N-1;
frequencies = [0.25, 0.5, 0.75]; % Normalized frequencies
amplitudes = [1, 2, 1.5];
phases = [0, pi/4, pi/2];
% Varying SNR
snrValues = [5, 10, 15, 20]; % Different SNR values to test
% Prepare a figure for the combined plot
figure;
hold on; % Enable plotting on the same axes
for i = 1:length(snrValues)
SNR = snrValues(i);
signal = amplitudes(1)*exp(1i*2*pi*frequencies(1)*n + phases(1)) + ...
amplitudes(2)*exp(1i*2*pi*frequencies(2)*n + phases(2)) + ...
amplitudes(3)*exp(1i*2*pi*frequencies(3)*n + phases(3));
signal = awgn(signal, SNR, 'measured'); % Add white Gaussian noise at specified SNR
% MUSIC algorithm
R = corrmtx(signal, 3, 'cov');
[Pxx, F] = pmusic(R, 3, 512, 1);
% Plotting on the same axes
plot(F, 10*log10(Pxx), 'LineWidth', 1.5);
end
hold off; % Disable plotting on the same axes
title('Combined MUSIC Algorithm Spectral Density at Different SNR');
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Spectral Power Density (dB)');
legend('SNR = 5 dB', 'SNR = 10 dB', 'SNR = 15 dB', 'SNR = 20 dB');
grid on;