% Read two gunshot recordings, plot waveforms, spectrums, and cross
% correlation

%% initialization
fs = 100000;
nfft1 = 2048;
nfft2 = 2048;
title1 = sprintf('%.0f point fft' ,nfft1);
title2 = sprintf('%.0f point fft' ,nfft2);
file1loc = 'Z:\jtobin\gunshots\PoliceAcademy1\recording_043_0.wav';
file2loc = 'Z:\jtobin\gunshots\PoliceAcademy1\recording_044_0.wav';

%% read and resample data
[y1, ~] = audioread(file1loc);
[y2, ~] = audioread(file2loc);
y1 = resample(y1,1,10); % downsample from 1Msps to 100ksps
y2 = resample(y2,1,10);
% y1 = resample(y1,1,40); % downsample from 1Msps to 40ksps
% y2 = resample(y2,1,40); % downsample from 1Msps to 40ksps

%sound(y1,fs);

%% generate indices
t1_index = transpose(linspace(1,length(y1),length(y1)));
t2_index = transpose(linspace(1,length(y2),length(y2)));
f1_index = transpose(fs/2*linspace(0,1,nfft1/2+1));
f2_index = transpose(fs/2*linspace(0,1,nfft2/2+1));

yfft1 = fft(y1,nfft1);
yfft2 = fft(y2,nfft2);

%% plots
subplot(5,1,1);
plot(t1_index,y1);
title('Sample 1');
grid minor;
xlabel('Sample index');
ylabel('Amplitude');

subplot(5,1,2);
plot(t2_index,y2);
title('Sample 2');
grid minor;
xlabel('Sample index');
ylabel('Amplitude');

subplot(5,1,3);
% plot(f1_index, abs(yfft1(1: nfft1/2+1)));
plot(f1_index(1:floor(2*length(f1_index)/5)),abs(yfft1(1:floor(2*length(f1_index)/5)))); 
% (1:floor(2*length(f1_index)/5)) is 2/5 of fs/2 = 20k if fs = 100k
title(title1);
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,1,4);
% plot(f2_index, abs(yfft1(1: nfft2/2+1)));
plot(f2_index(1:floor(2*length(f2_index)/5)), abs(yfft2(1:floor(2*length(f2_index)/5))));%nfft2/2+1)));
title(title2);
grid minor;
xlabel('Frequency');
ylabel('Amplitude');

a = xcorr(y1, y2);
a_index = linspace(1,length(a),length(a));

subplot(5,1,5);
plot(a_index, a);
title('Cross correlation');
grid minor;


