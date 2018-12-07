m4AFilename = 'voice.m4a';
[Dual_Channel_voice,Fs] = audioread(m4AFilename);
Single_Channel_voice = transpose(Dual_Channel_voice(:,1));

% AM
t = 1:length(Single_Channel_voice);
mu = 0.9;
Ac = 3;
fc = 250;
AM = Ac .* (1 + mu.*Single_Channel_voice) .* cos(2*pi*fc*t);
noisy_AM = wgn(1, length(Single_Channel_voice), -35) + AM;
AM_hilbert_noisy = imag(hilbert(noisy_AM));
envelop_AM_noisy = sqrt(AM_hilbert_noisy.^2 + noisy_AM.^2)/10;
AM_Dualband = zeros(length(Single_Channel_voice),2);
AM_Dualband(:,1) = envelop_AM_noisy;
AM_Dualband(:,2) = envelop_AM_noisy;
AM_File = 'AM.m4a';
audiowrite(AM_File,AM_Dualband,Fs)

% DSB
DSB = Ac .* Single_Channel_voice .* cos(2*pi*fc*t);
[b,a] = butter(5,fc/(3000/2)); % Butterworth filter of order 5
noisy_DSB = wgn(1, length(Single_Channel_voice), -30) + DSB;
DSB_Demodulated = 2*filter(b,a,noisy_DSB .* cos(2*pi*fc*t))./Ac;
DSB_Dualband = zeros(length(Single_Channel_voice),2);
DSB_Dualband(:,1) = DSB_Demodulated;
DSB_Dualband(:,2) = DSB_Demodulated;
DSB_File = 'DSB.m4a';
audiowrite(DSB_File,DSB_Dualband,Fs)

% SSB
Single_Channel_voice_Hilbert = imag(hilbert(Single_Channel_voice));
SSB = 0.5 .* Ac .* (Single_Channel_voice .* cos(2*pi*fc*t) + Single_Channel_voice_Hilbert .* sin(2*pi*fc*t));
noisy_SSB = wgn(1, length(Single_Channel_voice), -35) + SSB;
SSB_Demodulated = 4*filter(b,a,noisy_SSB .* cos(2*pi*fc*t))./Ac;
SSB_Dualband = zeros(length(Single_Channel_voice),2);
SSB_Dualband(:,1) = SSB_Demodulated;
SSB_Dualband(:,2) = SSB_Demodulated;
SSB_File = 'SSB.m4a';
audiowrite(SSB_File,SSB_Dualband,Fs)

% PM
frequency_deviation = 100;
phase_deviation = pi/4;
PM_1 = Ac*cos(2*pi*fc*t + phase_deviation * Single_Channel_voice);
PM = pmmod(Single_Channel_voice,fc, 600, phase_deviation);
noisy_PM = wgn(1, length(Single_Channel_voice), -50) + PM;
PM_Demodulated = pmdemod(noisy_PM,fc,600,phase_deviation);
PM_Dualband = zeros(length(Single_Channel_voice),2);
PM_Dualband(:,1) = PM_Demodulated;
PM_Dualband(:,2) = PM_Demodulated;
PM_File = 'PM.m4a';
audiowrite(PM_File,PM_Dualband,Fs)

% FM
FM_1 = Ac*cos(2*pi*fc*t + 2*pi*frequency_deviation * Single_Channel_voice);
FM = fmmod(Single_Channel_voice,fc, 600, frequency_deviation);
noisy_FM = wgn(1, length(Single_Channel_voice), -70) + FM;
FM_Demodulated = fmdemod(noisy_PM,fc,600,frequency_deviation);
FM_Dualband = zeros(length(Single_Channel_voice),2);
FM_Dualband(:,1) = FM_Demodulated;
FM_Dualband(:,2) = FM_Demodulated;
FM_File = 'FM.m4a';
audiowrite(FM_File,FM_Dualband,Fs)



