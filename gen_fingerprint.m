% Read gunshot recording and save set range as a fingerprint

clear all

%% initialization

% samplelength should always be equal to or less than nfft to allow for
% proper amount of zero padding
nfft = 4096;
samplelength = 2048;
file = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\R_27.wav';

%% Read and resample data. 

[y, fs] = audioread(file);

% (:,1) is first channel, (:,2) is second
y1 = y(:,1);
y2 = y(:,2);

% Downsample
fdes = 100000; % desired sampling frequency
y1 = resample(y1,fdes,fs);
y2 = resample(y2,fdes,fs);

%% prepare name string

% save as (sample name)_(shock/muz)_(ch1/2)_sampleLength_nfft_fs
shot_name_ind = strfind(file, '_');
shot_name = file(shot_name_ind-1 : shot_name_ind +2);

if fdes >99999
    c = 3;
else
    c = 2;
end

fs_str = num2str(fdes);
shock1name = strcat('Z:\jtobin\gunshots\fingerprintLib\', shot_name, '_s1_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.mat');
shock2name = strcat('Z:\jtobin\gunshots\fingerprintLib\', shot_name, '_s2_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.mat');
muz1name = strcat('Z:\jtobin\gunshots\fingerprintLib\', shot_name, '_m1_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.mat');
muz2name = strcat('Z:\jtobin\gunshots\fingerprintLib\', shot_name, '_m2_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.mat');

%% look at shockwave and muzzle blasts

% Shock wave arrives before muzzle blast if supersonic projectile passes by mic.
% Reflected sounds are included in the following pieces.

shockstart1 = 46500;
shockstart2 = 46400;
muzstart1 = 54000;
muzstart2 = 53900;

y1shock = y1(shockstart1:(shockstart1+samplelength-1));
y1muz = y1(muzstart1:(muzstart1+samplelength-1));
y2shock = y2(shockstart2:(shockstart2+samplelength-1));
y2muz = y2(muzstart2:(muzstart2+samplelength-1));

%% save sample 1 shockwave as a fingerprint

% zero pad sample to at least twice its length
zeropad = transpose(linspace(0, 0, nfft - samplelength));

y1shockpad = cat(1, y1shock, zeropad);
y1muzpad = cat(1, y1muz, zeropad);
y2shockpad = cat(1, y2shock, zeropad);
y2muzpad = cat(1, y2muz, zeropad);

y1s_fp = conj(fft(y1shockpad, nfft));
y1m_fp = conj(fft(y1muzpad, nfft));
y2s_fp = conj(fft(y2shockpad, nfft));
y2m_fp = conj(fft(y2muzpad, nfft));

savefiles = 1;

if savefiles == 0
    save(shock1name, 'y1s_fp');
%     save(muz1name, 'y1m_fp');
    save(shock2name, 'y2s_fp');
%     save(muz2name, 'y2m_fp');
end

%% generate indices

index_t = transpose(linspace(1,length(y1),length(y1)));
index_f1 = transpose(fdes/2*linspace(0,1,nfft/2));
index_f2 = transpose(fdes/2*linspace(0,1,nfft/2));
index_shock1 = transpose(linspace(shockstart1, shockstart1+samplelength, samplelength));
index_shock2 = transpose(linspace(shockstart2, shockstart2+samplelength, samplelength));
index_muz1 = transpose(linspace(muzstart1, muzstart1+samplelength, samplelength));
index_muz2 = transpose(linspace(muzstart2, muzstart2+samplelength, samplelength));

%% plot
% use to set m and n of subplot
a = 3;
b = 2;
figure;

% ch1
subplot(a,b,1);
plot(index_t,y1);
title('\bf Ch1');
grid minor;
xlabel('Sample Index');
ylabel('Amplitude');

% ch2
subplot(a,b,2);
plot(index_t,y2);
title('\bf Ch2');
grid minor;
xlabel('Sample Index');
ylabel('Amplitude');


% y1 shock
subplot(a,b,3);
plot(index_shock1,y1shock);
title('\bf Ch 1 Shockwave');
grid minor;
xlabel('Sample Index');
ylabel('Amplitude');
xlim([min(index_shock1) max(index_shock1)]);

% y1 muz
subplot(a,b,5);
plot(index_muz1,y1muz);
title('\bf Ch 1 Muzzleblast');
grid minor;
xlabel('Sample Index');
ylabel('Amplitude');
xlim([min(index_muz1) max(index_muz1)]);

% y2 shock
subplot(a,b,4);
plot(index_shock2,y2shock);
title('\bf Ch 2 Shockwave');
grid minor;
xlabel('Sample Index');
ylabel('Amplitude');
xlim([min(index_shock2) max(index_shock2)]);

% y2 muz
subplot(a,b,6);
plot(index_muz2,y2muz);
title('\bf Ch 2 Muzzleblast');
grid minor;
xlabel('Sample Index');
ylabel('Amplitude');
xlim([min(index_muz2) max(index_muz2)]);

% % y1 shock fft
% subplot(a,b,4);
% plot(index_croppedfft(1:floor(2*length(index_croppedfft)/5)), abs(y1shockfft(1:floor(2*length(index_croppedfft)/5))));
% title('\bf Sample 1 Shockwave FFT');
% grid minor;
% xlabel('Frequency');
% ylabel('Amplitude');