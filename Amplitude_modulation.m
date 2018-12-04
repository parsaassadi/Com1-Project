% Defining message signals
fs = 500;
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

% Adding Noise to signals
var(x1c_AM)
noise1 = wgn(1,2*h*fs+1,0);
x1c_AM = x1c_AM + noise1;
plot(t,x1c_AM)