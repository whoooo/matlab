% Read two gunshot recordings, plot waveforms, spectrums, and cross
% correlation

%% read in gunshot description file and store data


%% initialization
fs = 100000;
nfft1 = 4096;
nfft2 = 4096;
nfftcrop = 4096;
title1 = sprintf('%.0f  point fft of sample 1' , nfft1);
title2 = sprintf('%.0f   point fft of sample 2' , nfft2);
file1loc = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\Walther PPQ\X_32.wav';
file2loc = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\Walther PPQ\X_33.wav';

%% Read and resample data. 

[y1, fs] = audioread(file1loc);
[y2, fs] = audioread(file2loc);

% (:,1) is first channel, (:,2) is second
y1 = y1(:,1);
y2 = y2(:,1);

% Downsample
fdes = 100000; % desired frequency
y1 = resample(y1,fdes,fs);
y2 = resample(y2,fdes,fs);

% Play audio
%sound(y1,fs);
%sound(y2,fs);

% Shock wave arrive before muzzle blast if projectile passes by mic.
% Reflected sounds are included in the following pieces.
% Each piece should be the same length
shockstart1 = 105000;
shockstart2 = 93200;
muzstart1 = 13680;
muzstart2 = 4000;
croplength = 2048;
y1shock = y1(shockstart1:(shockstart1+croplength));
y1muz = y1(muzstart1:(muzstart1+croplength));
y2shock = y2(shockstart2:(shockstart2+croplength));
y2muz = y2(muzstart2:(muzstart2+croplength));

% y1 = resample(y1,1,40); % downsample from 1Msps to 40ksps
% y2 = resample(y2,1,40); % downsample from 1Msps to 40ksps

% zeropad = transpose(linspace(0, 0, length(y1)));
% y1 = cat(1,y1,zeropad);
% y2 = cat(1,y2,zeropad);

%% Take FFTs
y1fft = fft(y1,nfft1);
y2fft = fft(y2,nfft2);
y1shockfft = fft(y1shock, nfftcrop);
y1muzfft = fft(y1muz, nfftcrop);
y2shockfft = fft(y2shock, nfftcrop);
y2muzfft = fft(y2muz, nfftcrop);

%% Take Cross Correlation and Cross Covariance
corr = xcorr(y1, y2);
corrshock = xcorr(y1shock, y2shock);
corrmuz = xcorr(y1muz, y2muz);

crosscov = xcov(y1,y2);
crosscovshock = xcov(y1shock, y2shock);
crosscovmuz = xcov(y1muz, y2muz);

%% generate indices
index_t1 = transpose((1000/fdes).*linspace(1,length(y1),length(y1)));
index_t2 = transpose((1000/fdes).*linspace(1,length(y2),length(y2)));
index_f1 = transpose(fdes/2*linspace(0,1,nfft1/2));
index_f2 = transpose(fdes/2*linspace(0,1,nfft2/2));
index_shock1 = transpose((1000/fdes).*linspace(shockstart1, shockstart1+croplength+1, croplength+1));
index_shock2 = transpose((1000/fdes).*linspace(shockstart2, shockstart2+croplength+1, croplength+1));
index_muz1 = transpose((1000/fdes).*linspace(muzstart1, muzstart1+croplength+1, croplength+1));
index_muz2 = transpose((1000/fdes).*linspace(muzstart2, muzstart2+croplength+1, croplength+1));
index_croppedfft = transpose(fdes/2*linspace(0,1,nfftcrop/2));
index_corr = linspace(1,length(corr),length(corr));
index_corr_short = linspace(1,length(crosscovshock),length(crosscovshock));

%% plot full waveform, muzzleblast, shockwave, and fft of sample 1

% use to set m and n of subplot
a = 3;
b = 2;
figure;

% y1
subplot(a,b,1);
plot(index_t1,y1);
title('\bf Sample 1');
grid minor;
xlabel('Time (mS)');
ylabel('Amplitude');

% y1 fft
subplot(a,b,2);
plot(index_f1(1:floor(2*length(index_f1)/5)),abs(y1fft(1:floor(2*length(index_f1)/5)))); 
%(1:floor(2*length(f1_index)/5)) is 2/5 of fs/2 = 20k if fs = 100k
title(title1, 'FontWeight', 'bold');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

% y1 shock
subplot(a,b,3);
plot(index_shock1,y1shock);
title('\bf Sample 1 Shockwave');
grid minor;
xlabel('Time (mS)');
ylabel('Amplitude');
xlim([min(index_shock1) max(index_shock1)]);

% y1 shock fft
subplot(a,b,4);
plot(index_croppedfft(1:floor(2*length(index_croppedfft)/5)), abs(y1shockfft(1:floor(2*length(index_croppedfft)/5))));
title('\bf Sample 1 Shockwave FFT');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

% y1 muzzleblast
subplot(a,b,5);
plot(index_muz1,y1muz);
title('\bf Sample 1 Muzzleblast');
grid minor;
xlabel('Time (mS)');
ylabel('Amplitude');
xlim([min(index_muz1) max(index_muz1)]);

% y1 muzzleblast fft
subplot(a,b,6);
plot(index_croppedfft(1:floor(2*length(index_croppedfft)/5)), abs(y1muzfft(1:floor(2*length(index_croppedfft)/5))));
title('\bf Sample 1 Muzzleblast FFT');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

%% plot full waveform, muzzleblast, shockwave, and fft of sample 2

% use to set m and n of subplot
a = 3;
b = 2;
figure;

% y2
subplot(a,b,1);
plot(index_t2,y2);
title('\bf Sample 2');
grid minor;
xlabel('Time (mS)');
ylabel('Amplitude');

% y2 fft
subplot(a,b,2);
plot(index_f2(1:floor(2*length(index_f2)/5)), abs(y2fft(1:floor(2*length(index_f2)/5))));
%(1:floor(2*length(f1_index)/5)) is 2/5 of fs/2 = 20k if fs = 100k
title(title2, 'FontWeight', 'bold');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

% y2 shock
subplot(a,b,3);
plot(index_shock2,y2shock);
title('\bf Sample 2 Shockwave');
grid minor;
xlabel('Time (mS)');
ylabel('Amplitude');
xlim([min(index_shock2) max(index_shock2)]);

% y2 shock fft
subplot(a,b,4);
plot(index_croppedfft(1:floor(2*length(index_croppedfft)/5)), abs(y2shockfft(1:floor(2*length(index_croppedfft)/5))));
title('\bf Sample 2 Shockwave FFT');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

% y2 muzzleblast
subplot(a,b,5);
plot(index_muz2,y2muz);
title('\bf Sample 2 Muzzleblast');
grid minor;
xlabel('Time (mS)');
ylabel('Amplitude');
xlim([min(index_muz2) max(index_muz2)]);

% y2 muzzleblast fft
subplot(a,b,6);
plot(index_croppedfft(1:floor(2*length(index_croppedfft)/5)), abs(y2muzfft(1:floor(2*length(index_croppedfft)/5))));
title('\bf Sample 2 Muzzleblast FFT');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

%% plot cross correlation and cross covariance of complete waveforms, muzzleblast, and shockwave

figure

% use to set m and n of subplot
c = 3;
d = 1;

subplot(c,d,1);
plot(index_corr, corr);
title('\bf Cross Correlation of Full Signals');
grid minor;

% subplot(c,d,2);
% plot(index_corr, crosscov);
% title('\bf Cross Covariance of Full Signals');
% grid minor;

subplot(c,d,2);
plot(index_corr_short, corrshock);
title('\bf Cross Correlation of Shockwaves');
grid minor;

% subplot(c,d,4);
% plot(index_corr_short, crosscovshock);
% title('\bf Cross Covariance of Shockwaves');
% grid minor;

subplot(c,d,3);
plot(index_corr_short, corrmuz);
title('\bf Cross Correlation of Muzzleblast');
grid minor;

% subplot(c,d,6);
% plot(index_corr_short, crosscovmuz);
% title('\bf Cross Covariance of Muzzleblast');
% grid minor;

% subplot(6,2,);
% % plot(f1_index, abs(yfft1(1: nfft1/2+1)));
% plot(f1_index(1:floor(2*length(f1_index)/5)),abs(y1fft(1:floor(2*length(f1_index)/5)))); 
% % (1:floor(2*length(f1_index)/5)) is 2/5 of fs/2 = 20k if fs = 100k
% title(title1);
% grid minor;
% xlabel('Frequency');
% ylabel('Amplitude');
% 
% subplot(6,1,4);
% % plot(f2_index, abs(yfft1(1: nfft2/2+1)));
% plot(f2_index(1:floor(2*length(f2_index)/5)), abs(y2fft(1:floor(2*length(f2_index)/5))));%nfft2/2+1)));
% title(title2);
% grid minor;
% xlabel('Frequency');
% ylabel('Amplitude');
% 
% 
% 
% subplot(6,1,5);
% plot(a_index, corr);
% title('Cross correlation');
% grid minor;


