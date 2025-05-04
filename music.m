% Call the function with example parameters
pmusic_spectrum(200, [0.257, 0.2], [0, 0], 0.01, 4, 1024);
function pmusic_spectrum(signal_length, frequencies, phases, noise_variance,model_order, fft_points)
% Parameters
n = 0:signal_length-1; % Time vector
x = cos(frequencies(1)*pi*n + phases(1)) + sin(frequencies(2)*pi*n + phases(2)) + sqrt(noise_variance)*randn(size(n));
% Estimate the pseudospectrum using MUSIC
[pxx, fxx] = pmusic(x, model_order, fft_points, 1);
% Find peaks in the pseudospectrum
[pks, locs] = findpeaks(10*log10(pxx));
peak_freqs = fxx(locs);
% Plotting the pseudospectrum
figure;
plot(fxx, 10*log10(pxx), '-')
grid on
xlabel('Frequency (normalized)')
ylabel('Pseudospectrum (dB)')
title(sprintf('Pseudospectrum Estimate Using MUSIC (Model Order = %d)',model_order))
legend('Pseudospectrum')
end