% for use with xcorr_imp_test_loop_mult_fp (top, not topv2)

clear all;

% matlab automatically pads to next power of two size when doing corr
% scaling_factor = 2^sum(scaling_sch);
% fft_mat = fft(input);
% mag_mat = abs(fft_mat)/scaling_factor;

%% set n_fft, run/rst, ram initialized, and sampling frequency 
n_samples = 2048;
n_fft = n_samples * 2;
run = '1';
run_once = '0';
rst = '0';
use_adc = '1';
fs = 48000;
plot_result = 1;
plot_demo = 0;
plot_samp_and_fp = 1;
use_audiofile = 0;
use_fingerprint = 1;

%% get matlab generated results to compare (from audio file)

if use_audiofile == 1;
    
    sample = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\R_27.wav'; %s11
    [y, fs_orig] = audioread(sample);

    % (:,1) is first channel, (:,2) is second
    ych1 = y(:,1);
    ych2 = y(:,2);

    % Downsample
    ych1 = resample(ych1,fs,fs_orig);
    ych2 = resample(ych2,fs,fs_orig);

    % cut samples down
    shockstart1 = 46500;
    shockstart2 = 46400;
    y1shock = ych1(shockstart1:(shockstart1+n_samples-1));
    y2shock = ych2(shockstart2:(shockstart2+n_samples-1));

    % quantize and normalize time domain signal
    y1shock_quant = double(int16(y1shock/max(abs(y1shock)).*32768));
    y2shock_quant = double(int16(y2shock/max(abs(y2shock)).*32768));

    % zero pad sample to at least twice its length
    zeropad = transpose(linspace(0, 0, n_fft - n_samples));
    y1shockpad_quant = cat(1,y1shock_quant,zeropad);
    y2shockpad_quant = cat(1,y2shock_quant,zeropad);
    y1shockpad = cat(1,y1shock,zeropad);
    y2shockpad = cat(1,y2shock,zeropad);

    % take complex conjugate of fft
    y1shock_fp = conj(fft(y1shockpad, n_fft));
    y2shock_fp = conj(fft(y2shockpad, n_fft));

    % quantize the complex conjugate of the fft of quantized time domain signal
    y1shock_fp_quant = conj(fft(y1shockpad_quant, n_fft));
    y2shock_fp_quant = conj(fft(y2shockpad_quant, n_fft));
    y1shock_fp_quant = double(int16(y1shock_fp_quant/max(abs(y1shock_fp_quant)).*32768));
    y2shock_fp_quant = double(int16(y2shock_fp_quant/max(abs(y2shock_fp_quant)).*32768));

    % take manual cross correlation
    scaling_factor = 1;
    corr_man = fftshift(ifft((fft(y1shockpad, n_fft)./scaling_factor).*y2shock_fp));
    corr_man_quant = fftshift(ifft((fft(y1shockpad_quant, n_fft)./scaling_factor).*y2shock_fp_quant));
    corr_man_quant = double(int32(corr_man_quant/max(abs(corr_man_quant)).*2147483648));

    % use built in xcorr function
    corr_mat = xcorr(y1shock, y2shock);
    corr_mat_quant = xcorr(y1shock_quant, y2shock_quant);
    corr_mat_quant = double(int32(corr_mat_quant/max(abs(corr_mat_quant)).*2147483648));



end

%% get matlab generated results to compare (from fingerprint and time data files)

if use_fingerprint == 1
    
    % original time domain sample of fingerprint to use for matlab correlation
    fingerprint_t = transpose(load('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\A_33_s1_2048_48k.txt'));
    
    % fingerprint. conj(fft(fingerprint_t)) 
    fingerprint = load('Z:\jtobin\gunshots\fingerprintLib\f_domain\mat_files\A_33_s1_2048_4096_48k.txt'); % generated fingerprint file
    fingerprint = fingerprint(:,1) + 1j.*fingerprint(:,2); % combine real and imaginary components 
    
    % time domain signal which is combined with fingerprint for xcorr
    sample = transpose(load('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\R_27_s1_2048_48k.txt')); 
    sample_f = fft(sample, n_fft);                                                                           

    corr_man = ifftshift(ifft(sample_f.*fingerprint));
    corr_man = corr_man(2:n_fft);
    corr_mat = xcorr(sample, fingerprint_t);
             
end

%% generate indices

index_corr = linspace(1,n_fft-1,n_fft-1);
index_corr_mat = linspace(1,n_fft-1,n_fft-1);

%% convert values to binary representation

if (run == '1') && (rst == '1')
    run = '0';
    rst = '0';
    disp('Run and rst both set to 1');
end

% see page 45 of fft doc for values
if n_samples == 256
    n_points_b = '000';
elseif n_samples == 512
    n_points_b = '001';
elseif n_samples == 1024
    n_points_b = '010';
elseif n_samples == 2048
    n_points_b = '011';
elseif n_samples == 4096
    n_points_b = '100';
else
    n_points_b = '100';
end

cmd_config = bin2dec(strcat('0', '0', use_adc, run_once, n_points_b, '0'));
cmd = bin2dec(strcat(run, rst, use_adc, run_once, n_points_b, '0'));

%% serial comms
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (8 * n_fft)); % 64 bits x n_fft points
set(s, 'OutputBufferSize', 1);
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');
fopen(s);
fwrite(s, cmd_config, 'uint8');
fwrite(s, cmd, 'uint8');
sdata = fread(s, (n_fft*2), 'int32'); % 32 bits of [n_fft_points * 2 (real+imag components)] values
% cmd = bin2dec(strcat('0', ram_initialized, n_points_b, n_fft_b));
% fwrite(s,cmd,'uint8');
fclose(s);

%% split real/imag parts, take magnitude
a = 1;
b = 1;
sdata_im = zeros(1,n_fft);
sdata_re = zeros(1,n_fft);

for i = 1:2:(n_fft*2) 
    sdata_im(a) = sdata(i).*1j;
    a = a + 1;
end

for i = 2:2:(n_fft*2) 
    sdata_re(b) = sdata(i);
    b = b + 1;
end

corr_fpga = (sdata_re + sdata_im);
corr_fpga = corr_fpga(2:n_fft);

%% plot quantized and non quantized xcorrs

if plot_demo == 1
    figure;
    subplot(4,1,1);
    plot(index_corr, corr_man);
    title('Manual xcorr');
    subplot(4,1,2);
    plot(index_corr_mat, corr_mat);
    title('Built in xcorr');
    subplot(4,1,3);
    plot(index_corr, corr_man_quant);
    title('Manual xcorr- quantized');
    subplot(4,1,4);
    plot(index_corr_mat, corr_mat_quant);
    title('Built in xcorr- quantized');
end

%% plot result

% figure;
% subplot(3,1,1);
% plot(index_corr_mat, corr_mat_quant/scale);
% title('Built in xcorr- quantized');
% grid on;
% subplot(3,1,2);
% plot(index_corr, corr_man_quant/scale);
% title('Manual xcorr- quantized');
% grid on;
% subplot(3,1,3);
% plot(index_corr, flip(fftshift(corr_fpga)));

if plot_result == 1
    
    scale = 1;

    figure;
    subplot(3,1,1);
    plot(index_corr_mat, corr_mat./scale);
    title('Built in xcorr- quantized');
    grid on;
    subplot(3,1,2);
    plot(index_corr, corr_man./scale);
    title('Manual xcorr- quantized');
    grid on;
    subplot(3,1,3);
    plot(index_corr, fftshift(corr_fpga));
    grid on;
    title('FPGA xcorr');
    
    
    norm_corr_man = normalize(corr_man,1);
    norm_corr_fpga = fftshift(normalize(corr_fpga,1));  
    
    
    figure;
    plot(index_corr, norm_corr_man, 'b',...
        index_corr, norm_corr_fpga, 'r');
    legend('Matlab Correlation', 'FPGA Correlation');
    xlabel('Correlation Index');
    ylabel('Amplitude');
    grid on;
    
    abs_norm_corr_man = abs(normalize(corr_man,1));
    abs_norm_corr_fpga = abs(fftshift(normalize(corr_fpga,1)));
    
    
    figure;
    plot(index_corr, abs_norm_corr_man, 'b',...
        index_corr, abs_norm_corr_fpga, 'r');
    legend('Matlab Correlation', 'FPGA Correlation');
    xlabel('Correlation Index');
    ylabel('Amplitude');
    grid on;
    
    rmse_abs =sqrt(sum((abs_norm_corr_man(:)-abs_norm_corr_fpga(:)).^2)/numel(abs_norm_corr_man));
    rmse =sqrt(sum((norm_corr_man(:)-norm_corr_fpga(:)).^2)/numel(norm_corr_man));

    
end

