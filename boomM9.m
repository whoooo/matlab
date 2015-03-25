% load and analyze beretta m9 audio clip from boom library

%% initialization
sampleloc = 'Z:\jtobin\gunshots\berettam9.mp3';
nfft = 4096;
fs = 96000;

%% read file and take fft
% beretta9.mp3 has left and right audio channels
[zlong, fs] = audioread(sampleloc);
zlongch1 = zlong(:,1);
zlongch2 = zlong(:,2);
z = zlongch1(1:68000);
yfft = fft(z,nfft);

%% create indices
f_index = transpose(fs/2*linspace(0,1,nfft/2+1));
t_index = transpose(linspace(1,length(z),length(z)));
tlong_index = transpose(linspace(1,length(zlong),length(zlong)));

 sound(z,fs);

%% plot
% time domain
subplot(4,1,1);
plot(tlong_index,zlongch1);
title('\bf Beretta M9 Multiple Shots Ch1');
grid minor;
xlabel('Sample index');
ylabel('Amplitude');

% time domain
subplot(4,1,2);
plot(tlong_index,zlongch2);
title('\bf Beretta M9 Multiple Shots Ch2');
grid minor;
xlabel('Sample index');
ylabel('Amplitude');

% time domain
subplot(4,1,3);
plot(t_index,z);
title('\bf Beretta M9 Single Shot');
grid minor;
xlabel('Sample index');
ylabel('Amplitude');

% freq domain
subplot(4,1,4);
plot(f_index, 2*abs(yfft(1:nfft/2+1)));
title('\bf FFT of Single Shot');
grid minor;
xlabel('Frequency');
ylabel('Amplitude');
