# 📡 Subspace Frequency Estimation

This project investigates high-resolution subspace-based frequency estimation techniques using the **MUSIC (Multiple Signal Classification)** and **Eigenvector (EV)** methods. These algorithms are applied to identify frequencies embedded in noisy signals, comparing their effectiveness across different SNRs, amplitudes, and frequency spacings.

## 📘 Abstract

The study performs a comparative analysis of MUSIC and EV algorithms for spectral estimation. It simulates signals composed of complex exponentials with additive white Gaussian noise and evaluates algorithm performance using metrics like **Mean Squared Error (MSE)** and **Spectral Density** under various conditions.

## 🧠 Techniques Used

- **MUSIC Algorithm**  
- **Eigenvector Method**
- Signal modeling and noise injection
- Spectral analysis and MSE computation
- MATLAB for simulation and visualization

## 📊 Features

- Comparative performance analysis for MUSIC vs EV
- Variation of parameters: SNR, frequency spacing, amplitude levels
- Simulation of signal environments with both close and spread frequencies
- MSE plots and spectral density visualization

## 🔧 How to Run

1. Open the MATLAB script files in the `/code` or root folder.
2. Run simulations like:
   - `pmusic_spectrum(...)` to see MUSIC behavior
   - `peig(...)` to visualize the EV spectrum
3. Use included SNR and frequency variation loops to replicate comparison graphs

## 🖼️ Outputs

- Spectral density plots (MUSIC and EV)
- MSE vs. SNR plots for both algorithms
- Performance trends under different simulation conditions



