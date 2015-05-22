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

% save f domain as (sample name)_(shock/muz)(ch1/2)_sampleLength_nfft_fs
% save t domain as (sample name)_(shock/muz)(ch1/2)_sampleLength_fs

shot_name_ind = strfind(file, '_');
shot_name = file(shot_name_ind-1 : shot_name_ind +2);

if fdes >99999
    c = 3;
else
    c = 2;
end

fs_str = num2str(fdes);

shock1nameF = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\', shot_name, '_s1_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.coe');
shock2nameF = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\', shot_name, '_s2_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.coe');
muz1nameF = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\', shot_name, '_m1_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.coe');
muz2nameF = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\', shot_name, '_m2_', num2str(samplelength), '_', num2str(nfft), '_', fs_str(1:c), 'k', '.coe');

shock1nameT = strcat('Z:\jtobin\gunshots\fingerprintLib\time_domain\', shot_name, '_s1_', num2str(samplelength), '_', fs_str(1:c), 'k', '.coe');
shock2nameT = strcat('Z:\jtobin\gunshots\fingerprintLib\time_domain\', shot_name, '_s2_', num2str(samplelength), '_', fs_str(1:c), 'k', '.coe');

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

%% save fingerprint

zeropad = transpose(linspace(0, 0, nfft - samplelength));

% zero pad sample to at least twice its length
y1shockpad = cat(1, y1shock, zeropad);
y1muzpad = cat(1, y1muz, zeropad);
y2shockpad = cat(1, y2shock, zeropad);
y2muzpad = cat(1, y2muz, zeropad);

% take complex conjugate of fft
y1s_fp_f = conj(fft(y1shockpad, nfft));
y1m_fp_f = conj(fft(y1muzpad, nfft));
y2s_fp_f = conj(fft(y2shockpad, nfft));
y2m_fp_f = conj(fft(y2muzpad, nfft));

% normalize and quantize the complex conjugate of fft data
y1s_fp_norm = y1s_fp_f/abs(max(y1s_fp_f));
y1s_fp_quant = int16(y1s_fp_norm .*32768);
y2s_fp_norm = y2s_fp_f/abs(max(y2s_fp_f));
y2s_fp_quant = int16(y2s_fp_norm .*32768);

% concatenate real and imaginary parts of fingerprint
y1s_fp_real_bin = dec2bin(typecast(int16(real(y1s_fp_quant)),'uint16'));
y2s_fp_real_bin = dec2bin(typecast(int16(real(y2s_fp_quant)),'uint16'));
y1s_fp_imag_bin = dec2bin(typecast(int16(imag(y1s_fp_quant)),'uint16'));
y2s_fp_imag_bin = dec2bin(typecast(int16(imag(y2s_fp_quant)),'uint16'));

y1s_fp = cell2mat(strcat(y1s_fp_imag_bin, y1s_fp_real_bin, {' '}));
y2s_fp = cell2mat(strcat(y2s_fp_imag_bin, y2s_fp_real_bin, {' '}));
    
% %  compare normalized and quantized version of y1 fingerprint to original
% figure;
% subplot(2,1,1)
% plot(real(y1s_fp_f));
% subplot(2,1,2)
% plot(real(y1s_fp_quant));

% convert time domain sample to 16 bit int
y1s_t = int16(y1shock.*32768);
y2s_t = int16(y2shock.*32768);

savefiles = 1;

if savefiles == 1
    f_y1s_fp = fopen(shock1nameF, 'wt');
    fprintf(f_y1s_fp, '%s\n', 'memory_initialization_radix=2;', 'memory_initialization_vector= ');
    fprintf(f_y1s_fp, '%s ', transpose(y1s_fp));
    fprintf(f_y1s_fp, '%c', ';');
    fclose(f_y1s_fp);
    
    f_y2s_fp = fopen(shock2nameF, 'wt');
    fprintf(f_y2s_fp, '%s\n', 'memory_initialization_radix=2;', 'memory_initialization_vector= ');
    fprintf(f_y2s_fp, '%s\n', transpose(y2s_fp));
    fprintf(f_y2s_fp, '%c', ';');
    fclose(f_y2s_fp);
    
    f_y1s_t = fopen(shock1nameT, 'wt');
    fprintf(f_y1s_t, '%s\n', 'memory_initialization_radix=10;', 'memory_initialization_vector= ');
    fprintf(f_y1s_t, '%i ', y1s_t);
    fprintf(f_y1s_t, '%c', ';');
    fclose(f_y1s_t);
    
    f_y2s_t = fopen(shock2nameT, 'wt');
    fprintf(f_y2s_t, '%s\n', 'memory_initialization_radix=10;', 'memory_initialization_vector= ');
    fprintf(f_y2s_t, '%i ', y2s_t);
    fprintf(f_y2s_t, '%c', ';');
    fclose(f_y2s_t);   
end

%% save half of sample 1 spectrum as fingerprint

% % zero pad sample to at least twice its length
% zeropad = transpose(linspace(0, 0, nfft - samplelength));
% 
% y1shockpad = cat(1, y1shock, zeropad);
% 
% y1s_f = fft(y1shockpad, nfft);
% 
% y1s_f_half = y1s_f(1:2048);
% 
% y1s_fp = conj(y1s_f_half);
% 
% savefiles = 1;
% 
% if savefiles == 1
%     save('Z:\jtobin\gunshots\fingerprintLib\R_27_s1_1024_4096_100k_half_freq', 'y1s_fp');
% end


%% generate indices

index_t = transpose(linspace(1,length(y1),length(y1)));
index_f1 = transpose(fdes/2*linspace(0,1,nfft/2));
index_f2 = transpose(fdes/2*linspace(0,1,nfft/2));
index_shock1 = transpose(linspace(shockstart1, shockstart1+samplelength, samplelength));
index_shock2 = transpose(linspace(shockstart2, shockstart2+samplelength, samplelength));
index_muz1 = transpose(linspace(muzstart1, muzstart1+samplelength, samplelength));
index_muz2 = transpose(linspace(muzstart2, muzstart2+samplelength, samplelength));

%% plot scaled and unscaled ffts
% 
% y1_fft = fft(y1shock);
% y1_scaled_fft = fft(y1shock.*32768);
% figure;
% subplot(2,1,1)
% plot(index_f1(1:1024), y1_fft);
% subplot(2,1,2)
% plot(index_f1(1:1024), y1_scaled_fft);

% %% plot
% % use to set m and n of subplot
% a = 3;
% b = 2;
% figure;
% 
% % ch1
% subplot(a,b,1);
% plot(index_t,y1);
% title('\bf Ch1');
% grid minor;
% xlabel('Sample Index');
% ylabel('Amplitude');
% 
% % ch2
% subplot(a,b,2);
% plot(index_t,y2);
% title('\bf Ch2');
% grid minor;
% xlabel('Sample Index');
% ylabel('Amplitude');
% 
% 
% % y1 shock
% subplot(a,b,3);
% plot(index_shock1,y1shock);
% title('\bf Ch 1 Shockwave');
% grid minor;
% xlabel('Sample Index');
% ylabel('Amplitude');
% xlim([min(index_shock1) max(index_shock1)]);
% 
% % y1 muz
% subplot(a,b,5);
% plot(index_muz1,y1muz);
% title('\bf Ch 1 Muzzleblast');
% grid minor;
% xlabel('Sample Index');
% ylabel('Amplitude');
% xlim([min(index_muz1) max(index_muz1)]);
% 
% % y2 shock
% subplot(a,b,4);
% plot(index_shock2,y2shock);
% title('\bf Ch 2 Shockwave');
% grid minor;
% xlabel('Sample Index');
% ylabel('Amplitude');
% xlim([min(index_shock2) max(index_shock2)]);
% 
% % y2 muz
% subplot(a,b,6);
% plot(index_muz2,y2muz);
% title('\bf Ch 2 Muzzleblast');
% grid minor;
% xlabel('Sample Index');
% ylabel('Amplitude');
% xlim([min(index_muz2) max(index_muz2)]);

% % y1 shock fft
% subplot(a,b,4);
% plot(index_croppedfft(1:floor(2*length(index_croppedfft)/5)), abs(y1shockfft(1:floor(2*length(index_croppedfft)/5))));
% title('\bf Sample 1 Shockwave FFT');
% grid minor;
% xlabel('Frequency');
% ylabel('Amplitude');