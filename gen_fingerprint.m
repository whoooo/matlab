% Read gunshot recording and save set range as a fingerprint
clear all

%% initialization
% samplelength should always be equal to or less than n_fft to allow for
% proper amount of zero padding
n_fft = 4096;
samplelength = 2048;

save_custom = 0;
save_shots = 0;

plot_t_orig = 0; % original time domain signals
plot_custom = 1; % plot custom fingerprints
plot_spectrum = 0;
plot_quant_comparisons = 0;

file = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\Q_21.wav';

%% Read and resample data. 

[y, fs] = audioread(file);

% (:,1) is first channel, (:,2) is second
y1 = y(:,1);
y2 = y(:,2);

% Downsample
fdes = 48000; % desired sampling frequency
y1 = resample(y1,fdes,fs);
y2 = resample(y2,fdes,fs);

%% prepare name strings

% save f domain as (sample name)_(shock/muz)(ch1/2)_sampleLength_n_fft_fs
% save t domain as (sample name)_(shock/muz)(ch1/2)_sampleLength_fs

shot_name_ind = strfind(file, '_');
shot_name = file(shot_name_ind-1 : shot_name_ind +2);

if fdes >99999
    c = 3;
else
    c = 2;
end

fs_str = num2str(fdes);

shock1nameF = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\', shot_name, '_s1_', num2str(samplelength), '_', num2str(n_fft), '_', fs_str(1:c), 'k', '.coe');
shock2nameF = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\', shot_name, '_s2_', num2str(samplelength), '_', num2str(n_fft), '_', fs_str(1:c), 'k', '.coe');
shock1nameF_mat = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\mat_files\', shot_name, '_s1_', num2str(samplelength), '_', num2str(n_fft), '_', fs_str(1:c), 'k', '.txt');
shock2nameF_mat = strcat('Z:\jtobin\gunshots\fingerprintLib\f_domain\mat_files\', shot_name, '_s2_', num2str(samplelength), '_', num2str(n_fft), '_', fs_str(1:c), 'k', '.txt');

shock1nameT = strcat('Z:\jtobin\gunshots\fingerprintLib\time_domain\', shot_name, '_s1_', num2str(samplelength), '_', fs_str(1:c), 'k', '.coe');
shock2nameT = strcat('Z:\jtobin\gunshots\fingerprintLib\time_domain\', shot_name, '_s2_', num2str(samplelength), '_', fs_str(1:c), 'k', '.coe');
shock1nameT_mat = strcat('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\', shot_name, '_s1_', num2str(samplelength), '_', fs_str(1:c), 'k', '.txt');
shock2nameT_mat = strcat('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\', shot_name, '_s2_', num2str(samplelength), '_', fs_str(1:c), 'k', '.txt');

%% get shockwave and muzzle blasts

% Shock wave arrives before muzzle blast if supersonic projectile passes by mic.
% Reflected sounds are included in the following pieces.

shockstart1 = 48000;%46500*fdes/fs; 
shockstart2 = 48000;%*fdes/fs;

y1shock = y1(shockstart1:(shockstart1+samplelength-1));
y2shock = y2(shockstart2:(shockstart2+samplelength-1));

index_shock = linspace(0,length(y1shock),length(y1shock));

%% generate indices

index_t = transpose(linspace(1,length(y1),length(y1)));
index_f1 = transpose(fdes/2*linspace(0,1,n_fft/2));
index_f2 = transpose(fdes/2*linspace(0,1,n_fft/2));
index_fcomplete = transpose(fdes*linspace(0,1,n_fft));
index_shock1 = transpose(linspace(shockstart1, shockstart1+samplelength, samplelength));
index_shock2 = transpose(linspace(shockstart2, shockstart2+samplelength, samplelength));
% index_muz1 = transpose(linspace(muzstart1, muzstart1+samplelength, samplelength));
% index_muz2 = transpose(linspace(muzstart2, muzstart2+samplelength, samplelength));

%% plot original time domain signals

if plot_t_orig == 1 
    figure
    
    % plot entire sample with y1shock overlayed
    subplot(3,1,1)
    plot(index_t, y1, 'b', index_t(shockstart1:(shockstart1+samplelength)), y1(shockstart1:(shockstart1+samplelength)), 'r');
    title('Channel 1');
    
    % plot entire sample with y2shock overlayed
    subplot(3,1,2)
    plot(index_t, y2, 'b', index_t(shockstart2:(shockstart2+samplelength)), y2(shockstart2:(shockstart2+samplelength)), 'r');
    title('Channel 2');
    
   % plot shockwaves overlayed on each other
    subplot(3,1,3)
    plot(index_shock1, y1shock, 'b', index_shock2, y2shock, 'r');
    title('Channel 2');
    
end

%% create fingerprints

% normalize and quantize 
y1shock_quant = double(int16(y1shock/max(abs(y1shock)).*32767));
y2shock_quant = double(int16(y2shock/max(abs(y2shock)).*32767));

% zero pad sample to at least twice its length (up to n_fft in total length)
zeropad = transpose(linspace(0, 0, n_fft - samplelength));
y1shock_quant_pad = cat(1,y1shock_quant,zeropad);
y2shock_quant_pad = cat(1,y2shock_quant,zeropad);

% take complex conjugate of the fft of quantized time domain signal
y1shock_f = fft(y1shock_quant_pad, n_fft);
y2shock_f = fft(y2shock_quant_pad, n_fft);
y1shock_fp = conj(y1shock_f);
y2shock_fp = conj(y2shock_f);

% quantize f domain conjugate
y1shock_fp_quant = double(int16(y1shock_fp./max(abs(y1shock_fp)).*32767)); % add .
y2shock_fp_quant = double(int16(y2shock_fp./max(abs(y2shock_fp)).*32767));

% separate real and imaginary parts of fingerprint and convert to binary
y1shock_fp_real_bin = dec2bin(typecast(int16(real(y1shock_fp_quant)),'uint16'));
y2shock_fp_real_bin = dec2bin(typecast(int16(real(y2shock_fp_quant)),'uint16'));
y1shock_fp_imag_bin = dec2bin(typecast(int16(imag(y1shock_fp_quant)),'uint16'));
y2shock_fp_imag_bin = dec2bin(typecast(int16(imag(y2shock_fp_quant)),'uint16'));

% concatenate real and imaginary parts of fingerprint
y1s_fp = cell2mat(strcat(y1shock_fp_imag_bin, y1shock_fp_real_bin, {' '}));
y2s_fp = cell2mat(strcat(y2shock_fp_imag_bin, y2shock_fp_real_bin, {' '}));

%% plot spectrums

comp_y1_fft = abs(normalize(fft(y1shock, n_fft),1));
comp_y1_quant_fft = abs(normalize(fft(y1shock_quant, n_fft),1));
comp_y2_fft = abs(normalize(fft(y2shock, n_fft),1));
comp_y2_quant_fft = abs(normalize(fft(y2shock_quant, n_fft),1));

% rmse=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));;
rmse_ch1_fft =sqrt(sum((comp_y1_quant_fft(:)-comp_y1_fft(:)).^2)/numel(comp_y1_quant_fft));
rmse_ch2_fft =sqrt(sum((comp_y2_quant_fft(:)-comp_y2_fft(:)).^2)/numel(comp_y2_quant_fft));


if plot_spectrum == 1;
    
    figure
    
%     subplot(2,1,1)
    plot(index_fcomplete,comp_y1_fft, 'r');
    hold on;
    plot(index_fcomplete,comp_y1_quant_fft, 'b');
    xlim([0 20000]);
    legend('FFT of Original Signal', 'FFT of Quantized Signal');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    
    figure
    
%     subplot(2,1,2)
    plot(index_fcomplete, comp_y2_fft, 'r');
    hold on;
    plot(index_fcomplete,comp_y2_quant_fft, 'b');
    xlim([0 20000]);
    legend('FFT of Original Signal', 'FFT of Quantized Signal');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');

    figure
        
    subplot(4,1,1)
    plot(index_f1, abs(y1shock_f(1:length(y1shock_f)/2)));
    title('Y1 shock FFT');
    
    subplot(4,1,2)
    plot(index_f2, abs(y2shock_f(1:length(y2shock_f)/2)));
    title('Y2 shock FFT');
    
    subplot(4,1,3)
    plot(index_f1, abs(y1shock_fp_quant(1:length(y1shock_f)/2)));
    title('Y1 FP quant');
     
    subplot(4,1,4)
    plot(index_f2, abs(y2shock_fp_quant(1:length(y2shock_f)/2)));
    title('Y2 FP quant');   
        
    
    
end 
%% save fingerprints
    
if save_shots == 1
    
    % fp coe for y1
    f_y1s_fp = fopen(shock1nameF, 'wt');
    fprintf(f_y1s_fp, '%s\n', 'memory_initialization_radix=2;', 'memory_initialization_vector= ');
    fprintf(f_y1s_fp, '%s ', transpose(y1s_fp));
    fprintf(f_y1s_fp, '%c', ';');
    fclose(f_y1s_fp);
    
    % fp coe for y2
    f_y2s_fp = fopen(shock2nameF, 'wt');
    fprintf(f_y2s_fp, '%s\n', 'memory_initialization_radix=2;', 'memory_initialization_vector= ');
    fprintf(f_y2s_fp, '%s\n', transpose(y2s_fp));
    fprintf(f_y2s_fp, '%c', ';');
    fclose(f_y2s_fp);
    
    % t domain coe for y1
    f_y1s_t = fopen(shock1nameT, 'wt');
    fprintf(f_y1s_t, '%s\n', 'memory_initialization_radix=10;', 'memory_initialization_vector= ');
    fprintf(f_y1s_t, '%i ', y1shock_quant);
    fprintf(f_y1s_t, '%c', ';');
    fclose(f_y1s_t);
    
    % t domain coe for y2
    f_y2s_t = fopen(shock2nameT, 'wt');
    fprintf(f_y2s_t, '%s\n', 'memory_initialization_radix=10;', 'memory_initialization_vector= ');
    fprintf(f_y2s_t, '%i ', y2shock_quant);
    fprintf(f_y2s_t, '%c', ';');
    fclose(f_y2s_t);   
    
    % non binary y1 fingerprint file for use in matlab
%     y1shock_fp_quant_conj = conj(y1shock_fp_quant);
    f_y1s_fp_mat = fopen(shock1nameF_mat, 'wt');
    for i = 1:n_fft
        fprintf(f_y1s_fp_mat, '%e ', real(y1shock_fp_quant(i)));
        fprintf(f_y1s_fp_mat, '%e ', imag(y1shock_fp_quant(i)));  
        fprintf(f_y1s_fp_mat, '%s\n', ' ');
    end
    fclose(f_y1s_fp_mat);
    
    % non binary y2 fingerprint file for use in matlab
%     y2shock_fp_quant_conj = conj(y2shock_fp_quant);
    f_y2s_fp_mat = fopen(shock2nameF_mat, 'wt');
    for i = 1:n_fft
        fprintf(f_y2s_fp_mat, '%e ', real(y2shock_fp_quant(i)));
        fprintf(f_y2s_fp_mat, '%e ', imag(y2shock_fp_quant(i)));  
        fprintf(f_y2s_fp_mat, '%s\n', ' ');
    end
    fclose(f_y2s_fp_mat);
    
    % quantized non binary time domain y1
    f_y1s_t_mat = fopen(shock1nameT_mat, 'wt');
    fprintf(f_y1s_t_mat, '%i ', y1shock_quant);
    fclose(f_y1s_t_mat);
    
    % quantized non binary time domain y2
    f_y2s_t_mat = fopen(shock2nameT_mat, 'wt');
    fprintf(f_y2s_t_mat, '%i ', y2shock_quant);
    fclose(f_y2s_t_mat);       
end

%% save half length fingerprint (unfinished)

% % zero pad sample to at least twice its length
% zeropad = transpose(linspace(0, 0, n_fft - samplelength));
% 
% y1shockpad = cat(1, y1shock, zeropad);
% 
% y1s_f = fft(y1shockpad, n_fft);
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

%% create custom fingerprints (unfinished)
zeropad_cfp = linspace(0, 0, n_fft - samplelength);
index_custom = linspace(1,n_fft,n_fft);
index_custom_f = fdes/2*linspace(0,1,n_fft);

dirac_t = linspace(0,0,samplelength);
dirac_t(samplelength/2) = 65536;
dirac_t = cat(2,dirac_t,zeropad_cfp);
dirac_fft = fft(dirac_t, n_fft);
dirac_fp = conj(dirac_fft);
dirac_fp_quant = double(int16(dirac_fp./max(abs(dirac_fp)).*32768));


rect_t = linspace(0,0,samplelength);
rect_t(128:samplelength-128) = 65536;
rect_t = cat(2,rect_t,zeropad_cfp);
rect_fft = fft(rect_t, n_fft);
rect_fp = conj(rect_fft);
rect_fp_quant = double(int16(rect_fp./max(abs(rect_fp)).*32768));


if plot_custom == 1
    a = 4;
    b = 1;
    figure
    subplot(a,b,1)
    plot(index_custom, dirac_t);
    title('Dirac Delta');
    subplot(a,b,2)
    plot(index_custom_f(1:length(index_custom_f)/2), abs(dirac_fft(1:length(dirac_fft)/2)));
    title('FFT of Delta');
    subplot(a,b,3)
    plot(index_custom, rect_t);
    title('Boxcar');
    subplot(a,b,4)
    plot(index_custom_f(1:length(index_custom_f)/2), abs(rect_fft(1:length(rect_fft)/2)));  
    title('FFT of Boxcar');
end



% if save_custom == 1

%% comparison plots (unfinished)

if plot_quant_comparisons == 1
    
    y1shock_norm = normalize(y1shock,1);
    y1shock_quant_norm = normalize(y1shock_quant,1);
    y2shock_norm = normalize(y2shock,1);
    y2shock_quant_norm = normalize(y2shock_quant,1);   
    
    rmse_ch1_t =sqrt(sum((y1shock_norm(:)-y1shock_quant_norm(:)).^2)/numel(y1shock_norm));
    rmse_ch2_t =sqrt(sum((y2shock_norm(:)-y2shock_quant_norm(:)).^2)/numel(y2shock_norm));

      
    % compare quantized y1 signal to original (scaled to overlap)
%     scale_fact_t1 = max(abs(y1shock_quant))/max(abs(y1shock));
%     scaled_y1 = y1shock.* scale_fact_t1;
    figure
    plot(index_shock1, y1shock_norm, 'r', index_shock1, y1shock_quant_norm, 'b');
    legend('Original Signal', 'Quantized Signal')
    xlabel('Sample Index');
    ylabel('Amplitude');
%     title('Y1 T domain quantizattion comparison');
    
    % compare quantized y2 signal to original (scaled to overlap)
%     scale_fact_t2 = max(abs(y2shock_quant))/max(abs(y2shock));
%     scaled_y2 = y2shock.* scale_fact_t2;
    figure
    plot(index_shock2, y2shock_norm, 'r', index_shock2, y2shock_quant_norm, 'b');
    legend('Original Signal', 'Quantized Signal')
    xlabel('Sample Index');
    ylabel('Amplitude');
%     title('Y2 T domain quantizattion comparison');
    
    %  compare normalized and quantized version of y1 fingerprint to original
%     figure;
%     plot(index_f1(1:length(index_f1)/2),abs(y1shock_fp(1:length(y1shock_fp)/4)));
%     title('Y1S FFT');
    
    % compare fft of y1 to fingerprint 

% compare ifft of quantized fingerprint y1 to original quantized signal
    figure
    y1shock_fp_quant_ifft = ifft(y1shock_fp_quant, n_fft);
    scale_fact_ifft1 = max(abs(y1shock_fp_quant_ifft))/max(abs(y1shock_quant_pad));
    plot(linspace(0,n_fft,n_fft), flip(y1shock_fp_quant_ifft./scale_fact_ifft1), 'b', linspace(0,n_fft,n_fft), y1shock_quant_pad, 'r');
    title('IFFT of Y1 FP vs quantized Y1');
    legend('IFFT','Original');
    
% compare ifft of quantized fingerprint y2 to original quantized signal
    figure
    y2shock_fp_quant_ifft = ifft(y2shock_fp_quant, n_fft);
    scale_fact_ifft2 = max(abs(y2shock_fp_quant_ifft))/max(abs(y2shock_quant_pad));
    plot(linspace(0,n_fft,n_fft), flip(y2shock_fp_quant_ifft./scale_fact_ifft2), 'b', linspace(0,n_fft,n_fft), y2shock_quant_pad, 'r');
    title('IFFT of Y2 FP vs quantized Y2');
    legend('IFFT', 'Original');
    
end;