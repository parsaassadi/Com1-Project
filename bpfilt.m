function y = bpfilt(signal, f1, f2, fs, isplot)

if nargin < 4 || isempty(fs)
	fs = 1;
end

if nargin < 5 || isempty(isplot)
	isplot = 1;
end

if isrow(signal)
    signal = signal';
end
N  = length(signal);
dF = fs/N;
f  = (-fs/2:dF:fs/2-dF)';

if isempty(f1) || f1==-Inf
    BPF = (abs(f) < f2);
elseif isempty(f2) || f2==Inf
    BPF = (f1 < abs(f));
else
    BPF = ((f1 < abs(f)) & (abs(f) < f2));
end

%{
%% Power spectrum of the band-pass filter
if isplot
    figure;
    plot(f,BPF);
    title(sprintf('Power spectrum of the band-pass filter in (%.3f, %.3f) Hz',f1,f2));
end
%}

signal 	 = signal-mean(signal);
spektrum = fftshift(fft(signal))/N;

if isplot
    figure;
    subplot(2,1,1);
    plot(f,abs(spektrum));
    title('Power spectrum of the original signal');
end

spektrum = BPF.*spektrum;

if isplot
    subplot(2,1,2);
    plot(f,abs(spektrum));
    title(sprintf('Power spectrum of the band-pass filtered signal in (%.3f, %.3f) Hz',f1,f2));
end

y = ifft(ifftshift(spektrum)); %inverse ifft
y = real(y);

if isplot
    time = 1/fs*(0:N-1)';
	
    figure;
    subplot(2,1,1);
    plot(time,signal);
    title('The original time series');
    subplot(2,1,2);
    plot(time,y);
    title(sprintf('The band-pass filtered time series in (%.3f, %.3f) Hz',f1,f2));
end