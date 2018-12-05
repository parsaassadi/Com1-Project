% Defining message signals
fs = 5000;
h = 10;
t = -h:1/fs:h;
length(t)
x1 = (0).*(t<0) + (1).*(t>=0 & t<0.05) + (-2).*(t>=0.05 & t<0.1) + (0).*(t>=0.1);
x2 = (sinc(10*t).^2).*(t<=2 & t>=-2) + (0).*(t>2 & t<-2);


% AM modulation
mu = 0.4;
Ac1 = 3;
Ac2 = 1;
fc1 = 250;
fc2 = 100;
x1c_AM = Ac1 .* (1 + mu.*x1) .* cos(2*pi*fc1*t);
x2c_AM = Ac2 .* (1 + mu.*x2) .* cos(2*pi*fc2*t);

% DSB modulation
x1c_DSB = Ac1 .* x1 .* cos(2*pi*fc1*t);
x2c_DSB = Ac2 .* x2 .* cos(2*pi*fc2*t);

% SSB modulation
x1_hilbert = imag(hilbert(x1));
x2_hilbert = imag(hilbert(x2));
x1c_SSB = 0.5 .* Ac1 .* (x1 .* cos(2*pi*fc1*t) + x1_hilbert .* sin(2*pi*fc1*t));
x2c_SSB = 0.5 .* Ac2 .* (x2 .* cos(2*pi*fc2*t) + x2_hilbert .* sin(2*pi*fc2*t));

% Spectrum of message signals
f = fs/2*linspace(-1,1,2*h*fs+1);
fft_x1 = fftshift(fft(x1));
fft_x2 = fftshift(fft(x2));


% Spectrum of modulated signals
fft_x1c_AM = fftshift(fft(x1c_AM));
fft_x2c_AM = fftshift(fft(x2c_AM));
fft_x1c_DSB = fftshift(fft(x1c_DSB));
fft_x2c_DSB = fftshift(fft(x2c_DSB));
fft_x1c_SSB = fftshift(fft(x1c_SSB));
fft_x2c_SSB = fftshift(fft(x2c_SSB));

% DSB demodulation
[b,a] = butter(13,fc1/(fs/2)); % Butterworth filter of order 10
[c,d] = butter(10,fc2/(fs/2)); % Butterworth filter of order 10
x1_DSB_Demodulated = 2*filter(b,a,x1c_DSB .* cos(2*pi*fc1*t))./Ac1;
x2_DSB_Demodulated = 2*filter(c,d,x2c_DSB .* cos(2*pi*fc2*t))./Ac2;

% DSB demodulation - Asynchronous
phi0 = pi/6;
[b,a] = butter(13,fc1/(fs/2)); % Butterworth filter of order 10
[c,d] = butter(10,fc2/(fs/2)); % Butterworth filter of order 10
x1_DSB_Demodulated_Asynchronous = 2*filter(b,a,x1c_DSB .* cos(2*pi*fc1*t + phi0))./Ac1;
x2_DSB_Demodulated_Asynchronous = 2*filter(c,d,x2c_DSB .* cos(2*pi*fc2*t + phi0))./Ac2;

% SSB demodulation
x1_SSB_Demodulated = 4*filter(b,a,x1c_SSB .* cos(2*pi*fc1*t))./Ac1;
x2_SSB_Demodulated = 4*filter(c,d,x2c_SSB .* cos(2*pi*fc2*t))./Ac2;

% AM demodulation
x1c_AM_hilbert = imag(hilbert(x1c_AM));
envelope_x1c_AM = sqrt(x1c_AM_hilbert.^2 + x1c_AM.^2);
x2c_AM_hilbert = imag(hilbert(x2c_AM));
envelope_x2c_AM = sqrt(x2c_AM_hilbert.^2 + x2c_AM.^2);

% Noise addition and demodulating noisy signals
noisy_x1c_AM_1 = wgn(1, 2*h*fs+1, -25) + x1c_AM;
x1c_AM_hilbert_noisy_1 = imag(hilbert(noisy_x1c_AM_1));
envelope_x1c_AM_noisy_1 = sqrt(x1c_AM_hilbert_noisy_1.^2 + noisy_x1c_AM_1.^2);

noisy_x1c_AM_2 = wgn(1, 2*h*fs+1, -50) + x1c_AM;
x1c_AM_hilbert_noisy_2 = imag(hilbert(noisy_x1c_AM_2));
envelope_x1c_AM_noisy_2 = sqrt(x1c_AM_hilbert_noisy_2.^2 + noisy_x1c_AM_2.^2);

noisy_x2c_AM_1 = wgn(1, 2*h*fs+1, -25) + x2c_AM;
x2c_AM_hilbert_noisy_1 = imag(hilbert(noisy_x2c_AM_1));
envelope_x2c_AM_noisy_1 = sqrt(x2c_AM_hilbert_noisy_1.^2 + noisy_x2c_AM_1.^2);

noisy_x2c_AM_2 = wgn(1, 2*h*fs+1, -50) + x2c_AM;
x2c_AM_hilbert_noisy_2 = imag(hilbert(noisy_x2c_AM_2));
envelope_x2c_AM_noisy_2 = sqrt(x2c_AM_hilbert_noisy_2.^2 + noisy_x2c_AM_2.^2);

noisy_x1c_DSB_1 = wgn(1, 2*h*fs+1, -25) + x1c_DSB;
noisy_x1c_DSB_2 = wgn(1, 2*h*fs+1, -50) + x1c_DSB;
x1_DSB_Demodulated_noisy_1 = 2*filter(b,a,noisy_x1c_DSB_1 .* cos(2*pi*fc1*t))./Ac1;
x1_DSB_Demodulated_noisy_2 = 2*filter(b,a,noisy_x1c_DSB_2 .* cos(2*pi*fc1*t))./Ac1;


noisy_x2c_DSB_1 = wgn(1, 2*h*fs+1, -25) + x2c_DSB;
noisy_x2c_DSB_2 = wgn(1, 2*h*fs+1, -50) + x2c_DSB;
x2_DSB_Demodulated_noisy_1 = 2*filter(c,d,noisy_x2c_DSB_1 .* cos(2*pi*fc2*t))./Ac2;
x2_DSB_Demodulated_noisy_2 = 2*filter(c,d,noisy_x2c_DSB_2 .* cos(2*pi*fc2*t))./Ac2;

noisy_x1c_SSB_1 = wgn(1, 2*h*fs+1, -25) + x1c_SSB;
noisy_x1c_SSB_2 = wgn(1, 2*h*fs+1, -50) + x1c_SSB;
x1_SSB_Demodulated_noisy_1 = 4*filter(b,a,noisy_x1c_SSB_1 .* cos(2*pi*fc1*t))./Ac1;
x1_SSB_Demodulated_noisy_2 = 4*filter(b,a,noisy_x1c_SSB_2 .* cos(2*pi*fc1*t))./Ac1;


noisy_x2c_SSB_1 = wgn(1, 2*h*fs+1, -25) + x2c_SSB;
noisy_x2c_SSB_2 = wgn(1, 2*h*fs+1, -50) + x2c_SSB;
x2_SSB_Demodulated_noisy_1 = 4*filter(c,d,noisy_x2c_SSB_1 .* cos(2*pi*fc2*t))./Ac2;
x2_SSB_Demodulated_noisy_2 = 4*filter(c,d,noisy_x2c_SSB_2 .* cos(2*pi*fc2*t))./Ac2;

% Plotting section

