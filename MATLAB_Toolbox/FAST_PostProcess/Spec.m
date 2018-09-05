% Spectral analysis

% data = simout;
data = fastout;

% x = data.Wave1Elev;
x = data.TwrBsFxt;
% % x = data.TTDspFA;
% x = data.PtfmPitch; %(1001:end);
% x = data.GenSpeed;
% x = p2.data;

fs = 1/(data.Time(2) - data.Time(1));
% fs = 80;
N = length(x)-1;
xft = fft(x(2:end));
xft = xft(1:N/2+1);
psd = (1/(fs*N)) * abs(xft).^2;
psd(2:end-1) = 2*psd(2:end-1);
freq = 0:fs/N:fs/2;

% hold on
semilogx(freq,10*log10(psd))
hold on

%%

% Fs = 1000;
% t = 0:1/Fs:1-1/Fs;
% x = cos(2*pi*100*t); %+ randn(size(t));
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/N:Fs/2;
% 
% plot(freq,10*log10(psdx))
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')