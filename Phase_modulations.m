% Defining message
fs = 6000;
fm = 10;
h = 10;
t = -h:1/fs:h;
x = cos(20*pi*t).*(t<=2 & t>=0) + (0).*(t<0 & t>2);
integral_x = ((sin(20*pi*t))/(20*pi)).*(t<=2 & t>=0) + (0).*(t<0 & t>2);
frequency_deviation = 100;
phase_deviation = pi/4;
fc = 250;
Ac = 1;

% PM modulation
xc_pm = Ac*cos(2*pi*fc*t + phase_deviation * x);

% FM modulation
xc_fm = Ac*cos(2*pi*fc*t + 2*pi * frequency_deviation .* integral_x);

% Plotting message and modulated signals
figure
plot(t,x)
axis([-0.25 2.25 -2 2])
title('Message x')

figure
plot(t,xc_pm)
axis([-0.25 2.25 -2 2])
title('pm modulated without nosie')


figure
plot(t,xc_fm)
axis([-0.25 2.25 -2 2])
title('fm modulated without noise')

% Adding noise to modulated signals
noisy_xc_fm_10db = xc_fm + wgn(1,2*h*fs+1,-35);
noisy_xc_fm_20db = xc_fm + wgn(1,2*h*fs+1,-45);
noisy_xc_pm_10db = xc_pm + wgn(1,2*h*fs+1,-35);
noisy_xc_pm_20db = xc_pm + wgn(1,2*h*fs+1,-45);

figure
plot(t,noisy_xc_fm_10db)
axis([-0.25 2.25 -2 2])
title('fm modulated - snr:10db');

figure
plot(t,noisy_xc_fm_20db)
axis([-0.25 2.25 -2 2])
title('fm modulated - snr:20db');

figure
plot(t,noisy_xc_pm_10db)
axis([-0.25 2.25 -2 2])
title('pm modulated - snr:10db');

figure
plot(t,noisy_xc_pm_20db)
axis([-0.25 2.25 -2 2])
title('pm modulated - snr:20db');

% Spectrum analysis
f = fs/2*linspace(-1,1,2*h*fs+1);
fft_message = fftshift(fft(x));
fft_fm = fftshift(fft(xc_fm));
fft_pm = fftshift(fft(xc_pm));
fft_pm_10db = fftshift(fft(noisy_xc_pm_10db));
fft_pm_20db = fftshift(fft(noisy_xc_pm_20db));
fft_fm_10db = fftshift(fft(noisy_xc_fm_10db));
fft_fm_20db = fftshift(fft(noisy_xc_fm_20db));

figure
plot(f,abs(fft_message))
axis([-40 40 0 6000])
title('Message spectrum')

figure
plot(f,abs(fft_fm))
axis([-450 450 0 55000])
title('fm modulated spectrum - without noise')

figure
plot(f,abs(fft_fm_10db))
axis([-450 450 0 55000])
title('fm modulated spectrum - with 10db snr')

figure
plot(f,abs(fft_fm_20db))
axis([-450 450 0 55000])
title('fm modulated spectrum - with 20db snr')

figure
plot(f,abs(fft_pm))
axis([-450 450 0 60000])
title('pm modulated spectrum - without noise')

figure
plot(f,abs(fft_pm_10db))
axis([-450 450 0 60000])
title('pm modulated spectrum - with 10db snr')

figure
plot(f,abs(fft_pm_20db))
axis([-450 450 0 60000])
title('pm modulated spectrum - with 20db snr')

% Demodulating signals
fm_demod = fmdemod(xc_fm,fc,fs,frequency_deviation);
fm_demod_10db = fmdemod(noisy_xc_fm_10db,fc,fs,frequency_deviation);
fm_demod_20db = fmdemod(noisy_xc_fm_20db,fc,fs,frequency_deviation);

figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t,x)
xlim(ax1,[-0.25 2.25])
ylim(ax1,[-2 2])
title(ax1,'Message x')
plot(ax2,t,fm_demod)
xlim(ax2,[-0.25 2.25])
ylim(ax2,[-2 2])
title(ax2,'FM demodulated x without noise')

figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t,x)
xlim(ax1,[-0.25 2.25])
ylim(ax1,[-2 2])
title(ax1,'Message x')
plot(ax2,t,fm_demod_10db)
xlim(ax2,[-0.25 2.25])
ylim(ax2,[-2 2])
title(ax2,'FM demodulated x with 10db snr')

figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t,x)
xlim(ax1,[-0.25 2.25])
ylim(ax1,[-2 2])
title(ax1,'Message x')
plot(ax2,t,fm_demod_20db)
xlim(ax2,[-0.25 2.25])
ylim(ax2,[-2 2])
title(ax2,'FM demodulated x with 20db snr')


pm_demod = pmdemod(xc_pm,fc,fs,phase_deviation);
pm_demod_10db = pmdemod(noisy_xc_pm_10db,fc,fs,phase_deviation);
pm_demod_20db = pmdemod(noisy_xc_pm_20db,fc,fs,phase_deviation);

figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t,x)
xlim(ax1,[-0.25 2.25])
ylim(ax1,[-2 2])
title(ax1,'Message x')
plot(ax2,t,pm_demod)
xlim(ax2,[-0.25 2.25])
ylim(ax2,[-2 2])
title(ax2,'PM demodulated x without noise')

figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t,x)
xlim(ax1,[-0.25 2.25])
ylim(ax1,[-2 2])
title(ax1,'Message x')
plot(ax2,t,pm_demod_10db)
xlim(ax2,[-0.25 2.25])
ylim(ax2,[-2 2])
title(ax2,'PM demodulated x with 10db snr')

figure
ax1 = subplot(2,1,1); % top subplot
ax2 = subplot(2,1,2); % bottom subplot
plot(ax1,t,x)
xlim(ax1,[-0.25 2.25])
ylim(ax1,[-2 2])
title(ax1,'Message x')
plot(ax2,t,pm_demod_20db)
xlim(ax2,[-0.25 2.25])
ylim(ax2,[-2 2])
title(ax2,'PM demodulated x with 20db snr')



