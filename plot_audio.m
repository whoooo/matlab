clear all

%% get audio file

fs = 44100;
nfft = 32768;
sample_start = 1;

% [y, fs_orig] = audioread('Z:\jtobin\other_sounds\breaking-glass-1.wav');
% [y, fs_orig] = audioread('Z:\jtobin\other_sounds\breaking_window_indoors.wav');
[y, fs_orig] = audioread('Z:\jtobin\other_sounds\breaking_glass_sheet1.flac');

% plot_title = 'breaking-glass-1.wav';
% plot_title = 'breaking window indoors.wav';
plot_title = 'breaking glass sheet1.flac';

ych1 = resample(y(:,1), fs, fs_orig);
% ych2 = resample(y(:,2), fs, fs_orig);

ych1 = ych1(sample_start:sample_start + nfft);
% ych2 = ych2(sample_start:sample_start + nfft);

ych1_fft = fft(ych1, nfft);
% ych2_fft = fft(ych2, nfft);

index_t = 1000/fs.*linspace(1,length(ych1),length(ych1));
index_f = transpose(fs/2*linspace(0,1,nfft));

%% plot
figure
plot(index_t, ych1);
title(plot_title);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude');

set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

figure
plot(index_f(1:nfft/2), abs(normalize(ych1_fft(1:nfft/2),1)));
title(plot_title);
grid on;
xlabel('Frequency (Hz)');
ylabel('Amplitude');

figureHandle = gcf;

set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)


% figure
% plot(index_t, ych2);
% title('\bf Channel 2');
% grid on;
% xlabel('Time (ms)');
% ylabel('Amplitude');
% 
% figure
% plot(index_f(1:nfft/2), abs(ych2_fft(1:nfft/2)));
% title('\bf Channel 1');
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');