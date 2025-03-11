clc; clear; close all;
%Real-Time Signal Processing-Adaptive Noise Filtering & Fourier Analysis
% Parameters
fs = 1000; % Sampling frequency (Hz)
duration = 5; % Duration of the signal (seconds)
t = 0:1/fs:duration; % Time vector
f_signal = 50; % Signal frequency (Hz)
noise_level = 0.5;

% Initialize figure
figure;
h1 = subplot(3,1,1);
h2 = subplot(3,1,2);
h3 = subplot(3,1,3);

for k = 1:length(t)-fs
    % Generate a live signal with noise
    signal = sin(2*pi*f_signal*t(k:k+fs-1)) + noise_level * randn(1, fs);
    
    % Compute Fourier Transform
    N = length(signal);
    freq = (-N/2:N/2-1)*(fs/N); % Frequency vector
    fft_signal = fftshift(fft(signal));
    
    % Adaptive Filtering using LMS Algorithm
    mu = 0.01; % Learning rate
    order = 20; % Filter order
    b = zeros(order, 1); % Filter coefficients
    x = [zeros(order-1,1); signal']';
    for i = order:length(x)
        x_i = x(i:-1:i-order+1)';
        y = b' * x_i;
        e = signal(i-order+1) - y;
        b = b + mu * e * x_i;
    end
    filtered_signal = filter(b, 1, signal);
    
    % Update plots
    subplot(h1);
    plot(t(k:k+fs-1), signal);
    title('Live Noisy Signal'); xlabel('Time (s)'); ylabel('Amplitude');
    
    subplot(h2);
    plot(freq, abs(fft_signal)/N);
    title('Live Magnitude Spectrum'); xlabel('Frequency (Hz)'); ylabel('Magnitude');
    
    subplot(h3);
    plot(t(k:k+fs-1), filtered_signal);
    title('Live Filtered Signal'); xlabel('Time (s)'); ylabel('Amplitude');
    
    sgtitle('Real-Time Signal Processing: Fourier Transform & Adaptive Filtering');
    pause(0.1);
end
