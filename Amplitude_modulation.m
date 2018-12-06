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

% Noise addition and demodulating noisy signals and spectrum of noisy
% signals
noisy_x1c_AM_1 = wgn(1, 2*h*fs+1, -3.5) + x1c_AM;
x1c_AM_hilbert_noisy_1 = imag(hilbert(noisy_x1c_AM_1));
envelope_x1c_AM_noisy_1 = sqrt(x1c_AM_hilbert_noisy_1.^2 + noisy_x1c_AM_1.^2);

noisy_x1c_AM_2 = wgn(1, 2*h*fs+1, -13.46) + x1c_AM;
x1c_AM_hilbert_noisy_2 = imag(hilbert(noisy_x1c_AM_2));
envelope_x1c_AM_noisy_2 = sqrt(x1c_AM_hilbert_noisy_2.^2 + noisy_x1c_AM_2.^2);

noisy_x2c_AM_1 = wgn(1, 2*h*fs+1, -13) + x2c_AM;
x2c_AM_hilbert_noisy_1 = imag(hilbert(noisy_x2c_AM_1));
envelope_x2c_AM_noisy_1 = sqrt(x2c_AM_hilbert_noisy_1.^2 + noisy_x2c_AM_1.^2);

noisy_x2c_AM_2 = wgn(1, 2*h*fs+1, -23) + x2c_AM;
x2c_AM_hilbert_noisy_2 = imag(hilbert(noisy_x2c_AM_2));
envelope_x2c_AM_noisy_2 = sqrt(x2c_AM_hilbert_noisy_2.^2 + noisy_x2c_AM_2.^2);

noisy_x1c_DSB_1 = wgn(1, 2*h*fs+1, -22.5) + x1c_DSB;
noisy_x1c_DSB_2 = wgn(1, 2*h*fs+1, -32.5) + x1c_DSB;
x1_DSB_Demodulated_noisy_1 = 2*filter(b,a,noisy_x1c_DSB_1 .* cos(2*pi*fc1*t))./Ac1;
x1_DSB_Demodulated_noisy_2 = 2*filter(b,a,noisy_x1c_DSB_2 .* cos(2*pi*fc1*t))./Ac1;


noisy_x2c_DSB_1 = wgn(1, 2*h*fs+1, -37.7) + x2c_DSB;
noisy_x2c_DSB_2 = wgn(1, 2*h*fs+1, -47.6) + x2c_DSB;
x2_DSB_Demodulated_noisy_1 = 2*filter(c,d,noisy_x2c_DSB_1 .* cos(2*pi*fc2*t))./Ac2;
x2_DSB_Demodulated_noisy_2 = 2*filter(c,d,noisy_x2c_DSB_2 .* cos(2*pi*fc2*t))./Ac2;

noisy_x1c_SSB_1 = wgn(1, 2*h*fs+1, -26) + x1c_SSB;
noisy_x1c_SSB_2 = wgn(1, 2*h*fs+1, -36) + x1c_SSB;
x1_SSB_Demodulated_noisy_1 = 4*filter(b,a,noisy_x1c_SSB_1 .* cos(2*pi*fc1*t))./Ac1;
x1_SSB_Demodulated_noisy_2 = 4*filter(b,a,noisy_x1c_SSB_2 .* cos(2*pi*fc1*t))./Ac1;


noisy_x2c_SSB_1 = wgn(1, 2*h*fs+1, -41) + x2c_SSB;
noisy_x2c_SSB_2 = wgn(1, 2*h*fs+1, -50.7) + x2c_SSB;
x2_SSB_Demodulated_noisy_1 = 4*filter(c,d,noisy_x2c_SSB_1 .* cos(2*pi*fc2*t))./Ac2;
x2_SSB_Demodulated_noisy_2 = 4*filter(c,d,noisy_x2c_SSB_2 .* cos(2*pi*fc2*t))./Ac2;

fft_x1c_AM_noisy_1 = fftshift(fft(noisy_x1c_AM_1));
fft_x1c_AM_noisy_2 = fftshift(fft(noisy_x1c_AM_2));
fft_x2c_AM_noisy_1 = fftshift(fft(noisy_x2c_AM_1));
fft_x2c_AM_noisy_2 = fftshift(fft(noisy_x2c_AM_2));
fft_x1c_DSB_noisy_1 = fftshift(fft(noisy_x1c_DSB_1));
fft_x1c_DSB_noisy_2 = fftshift(fft(noisy_x1c_DSB_2));
fft_x2c_DSB_noisy_1 = fftshift(fft(noisy_x2c_DSB_1));
fft_x2c_DSB_noisy_2 = fftshift(fft(noisy_x2c_DSB_2));
fft_x1c_SSB_noisy_1 = fftshift(fft(noisy_x1c_SSB_1));
fft_x1c_SSB_noisy_2 = fftshift(fft(noisy_x1c_SSB_2));
fft_x2c_SSB_noisy_1 = fftshift(fft(noisy_x2c_SSB_1));
fft_x2c_SSB_noisy_2 = fftshift(fft(noisy_x2c_SSB_2));

% Plotting section
figure
plot(t, x1, 'LineWidth',2)
axis([-0.1 0.15 -3 1.5])
title('x1(t)')

figure
plot(t, x2,'LineWidth',2)
axis([-2.5 2.5 -1 1.5])
title('x2(t)')

figure
plot(t, x2c_AM,'LineWidth',2)
axis([-1.5 1.5 -1.5 1.5])
title('x2c(t)-Am modulated without noise')

figure
plot(t, x1c_AM,'LineWidth',2)
axis([-0.1 0.1 -7 7])
title('x1c(t)-Am modulated without noise')

figure
plot(t, x1c_DSB,'LineWidth',2)
axis([-0.15 0.15 -7 7])
title('x1c(t)-DSB modulated without noise')

figure
plot(t, x2c_DSB,'LineWidth',2)
axis([-0.5 0.5 -3 3])
title('x2c(t)-DSB modulated without noise')

figure
plot(t, x2c_SSB,'LineWidth',2)
axis([-0.5 0.5 -3 3])
title('x2c(t)-SSB modulated without noise')

figure
plot(t, x1c_SSB,'LineWidth',2)
axis([-0.15 0.3 -7 7])
title('x1c(t)-SSB modulated without noise')



figure
plot(t, noisy_x2c_AM_1,'LineWidth',2)
axis([-1.5 1.5 -1.5 1.5])
title('x2c(t)-Am modulated with 0.1 noise')

figure
plot(t, noisy_x1c_AM_1,'LineWidth',2)
axis([-0.1 0.1 -7 7])
title('x1c(t)-Am modulated with 0.1 noise')

figure
plot(t, noisy_x1c_DSB_1,'LineWidth',2)
axis([-0.15 0.15 -7 7])
title('x1c(t)-DSB modulated with 0.1 noise')

figure
plot(t, noisy_x2c_DSB_1,'LineWidth',2)
axis([-0.5 0.5 -3 3])
title('x2c(t)-DSB modulated with 0.1 noise')

figure
plot(t, noisy_x2c_SSB_1,'LineWidth',2)
axis([-0.5 0.5 -3 3])
title('x2c(t)-SSB modulated with 0.1 noise')

figure
plot(t, noisy_x1c_SSB_1,'LineWidth',2)
axis([-0.15 0.3 -7 7])
title('x1c(t)-SSB modulated with 0.1 noise')



figure
plot(t, noisy_x2c_AM_2,'LineWidth',2)
axis([-1.5 1.5 -1.5 1.5])
title('x2c(t)-Am modulated with 0.01 noise')

figure
plot(t, noisy_x1c_AM_2,'LineWidth',2)
axis([-0.1 0.1 -7 7])
title('x1c(t)-Am modulated with 0.01 noise')

figure
plot(t, noisy_x1c_DSB_2,'LineWidth',2)
axis([-0.15 0.15 -7 7])
title('x1c(t)-DSB modulated with 0.01 noise')

figure
plot(t, noisy_x2c_DSB_2,'LineWidth',2)
axis([-0.5 0.5 -3 3])
title('x2c(t)-DSB modulated with 0.01 noise')

figure
plot(t, noisy_x2c_SSB_2,'LineWidth',2)
axis([-0.5 0.5 -3 3])
title('x2c(t)-SSB modulated with 0.01 noise')

figure
plot(t, noisy_x1c_SSB_2,'LineWidth',2)
axis([-0.15 0.3 -7 7])
title('x1c(t)-SSB modulated with 0.01 noise')


figure
plot(f, abs(fft_x1), 'LineWidth',2)
axis([-400 400 0 700])
title('x1(t)-frequency spectrum')

figure
plot(f, abs(fft_x2),'LineWidth',2)
axis([-150 150 0 500])
title('x2(t)-frequency spectrum')

figure
plot(f, abs(fft_x1c_AM),'LineWidth',2)
axis([-260 260 0 160000])
title('x1c(t)-Am frequency spectrum without noise')

figure
plot(f, abs(fft_x2c_AM),'LineWidth',2)
axis([-120 120 0 60000])
title('x2c(t)-Am frequency spectrum without noise')

figure
plot(f, abs(fft_x1c_DSB),'LineWidth',2)
axis([-600 600 0 900])
title('x1c(t)-DSB frequency spectrum without noise')

figure
plot(f, abs(fft_x2c_DSB),'LineWidth',2)
axis([-200 200 0 400])
title('x2c(t)-DSB frequency spectrum without noise')

figure
plot(f, abs(fft_x2c_SSB),'LineWidth',2)
axis([-200 200 0 400])
title('x2c(t)-SSB frequency spectrum without noise')

figure
plot(f, abs(fft_x1c_SSB),'LineWidth',2)
axis([-500 500 0 900])
title('x1c(t)-SSB frequency spectrum without noise')




figure
plot(f, abs(fft_x1c_AM_noisy_1),'LineWidth',2)
axis([-260 260 0 160000])
title('x1c(t)-Am frequency spectrum with 0.1 noise')

figure
plot(f, abs(fft_x2c_AM_noisy_1),'LineWidth',2)
axis([-120 120 0 60000])
title('x2c(t)-Am frequency spectrum with 0.1 noise')

figure
plot(f, abs(fft_x1c_DSB_noisy_1),'LineWidth',2)
axis([-600 600 0 900])
title('x1c(t)-DSB frequency spectrum with 0.1 noise')

figure
plot(f, abs(fft_x2c_DSB_noisy_1),'LineWidth',2)
axis([-200 200 0 400])
title('x2c(t)-DSB frequency spectrum with 0.1 noise')

figure
plot(f, abs(fft_x2c_SSB_noisy_1),'LineWidth',2)
axis([-200 200 0 400])
title('x2c(t)-SSB frequency spectrum with 0.1 noise')

figure
plot(f, abs(fft_x1c_SSB_noisy_1),'LineWidth',2)
axis([-500 500 0 900])
title('x1c(t)-SSB frequency spectrum with 0.1 noise')
figure
plot(f, abs(fft_x1c_AM_noisy_2),'LineWidth',2)
axis([-260 260 0 160000])
title('x1c(t)-Am frequency spectrum with 0.01 noise')

figure
plot(f, abs(fft_x2c_AM_noisy_2),'LineWidth',2)
axis([-120 120 0 60000])
title('x2c(t)-Am frequency spectrum with 0.01 noise')

figure
plot(f, abs(fft_x1c_DSB_noisy_2),'LineWidth',2)
axis([-600 600 0 900])
title('x1c(t)-DSB frequency spectrum with 0.01 noise')

figure
plot(f, abs(fft_x2c_DSB_noisy_2),'LineWidth',2)
axis([-200 200 0 400])
title('x2c(t)-DSB frequency spectrum with 0.01 noise')

figure
plot(f, abs(fft_x2c_SSB_noisy_2),'LineWidth',2)
axis([-200 200 0 400])
title('x2c(t)-SSB frequency spectrum with 0.01 noise')

figure
plot(f, abs(fft_x1c_SSB_noisy_2),'LineWidth',2)
axis([-500 500 0 900])
title('x1c(t)-SSB frequency spectrum with 0.01 noise')
% 
figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,envelope_x1c_AM)
xlim(ax2,[-0.02 0.14])
title(ax2,'AM demodulated x1 without noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,envelope_x1c_AM_noisy_1)
xlim(ax2,[-0.02 0.14])
title(ax2,'AM demodulated x1 with 0.1 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,envelope_x1c_AM_noisy_2)
xlim(ax2,[-0.02 0.14])
title(ax2,'AM demodulated x1 with 0.01 noise')



figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t, x1_DSB_Demodulated)
xlim(ax2,[-0.02 0.14])
title(ax2,'DSB demodulated x1 without noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,x1_DSB_Demodulated_noisy_1)
xlim(ax2,[-0.02 0.14])
title(ax2,'DSB demodulated x1 with 0.1 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,x1_DSB_Demodulated_noisy_2)
xlim(ax2,[-0.02 0.14])
title(ax2,'DSB demodulated x1 with 0.01 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,x1_DSB_Demodulated_Asynchronous)
xlim(ax2,[-0.02 0.14])
title(ax2,'DSB demodulated x1 asynchronous (phi = pi/6)')



figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t, x1_SSB_Demodulated)
xlim(ax2,[-0.02 0.14])
title(ax2,'SSB demodulated x1 without noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,x1_SSB_Demodulated_noisy_1)
xlim(ax2,[-0.02 0.14])
title(ax2,'SSB demodulated x1 with 0.1 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x1);
xlim(ax1,[-0.02 0.14])
ylim(ax1,[-2.5 1.5])
title(ax1,'Message x1')
plot(ax2,t,x1_SSB_Demodulated_noisy_2)
xlim(ax2,[-0.02 0.14])
title(ax2,'SSB demodulated x1 with 0.01 noise')



figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,envelope_x2c_AM)
xlim(ax2,[-2.5 2.5])
title(ax2,'AM demodulated x2 without noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,envelope_x2c_AM_noisy_1)
xlim(ax2,[-2.5 2.5])
title(ax2,'AM demodulated x2 with 0.1 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,envelope_x2c_AM_noisy_2)
xlim(ax2,[-2.5 2.5])
title(ax2,'AM demodulated x2 with 0.01 noise')



figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t, x2_DSB_Demodulated)
xlim(ax2,[-2.5 2.5])
title(ax2,'DSB demodulated x2 without noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,x2_DSB_Demodulated_noisy_1)
xlim(ax2,[-2.5 2.5])
title(ax2,'DSB demodulated x2 with 0.1 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,x2_DSB_Demodulated_noisy_2)
xlim(ax2,[-2.5 2.5])
title(ax2,'DSB demodulated x2 with 0.01 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,x2_DSB_Demodulated_Asynchronous)
xlim(ax2,[-2.5 2.5])
title(ax2,'DSB demodulated x2 asynchronous (phi = pi/6)')



figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t, x2_SSB_Demodulated)
xlim(ax2,[-2.5 2.5])
title(ax2,'SSB demodulated x2 without noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,x2_SSB_Demodulated_noisy_1)
xlim(ax2,[-2.5 2.5])
title(ax2,'SSB demodulated x2 with 0.1 noise')

figure % new figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t, x2);
xlim(ax1,[-2.5 2.5])
ylim(ax1,[-0.1 1.5])
title(ax1,'Message x2')
plot(ax2,t,x2_SSB_Demodulated_noisy_2)
xlim(ax2,[-2.5 2.5])
title(ax2,'SSB demodulated x2 with 0.01 noise')
