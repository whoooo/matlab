clear all;

% matlab automatically pads to next power of two size when doing corr
% scaling_factor = 2^sum(scaling_sch);
% fft_mat = fft(input);
% mag_mat = abs(fft_mat)/scaling_factor;

%% set n_fft, run/rst, ram initialized, and sampling frequency 
n_samples = 2048;
n_fft = 4096;
run = '1'; % reset = 1 if run = 0
ram_initialized = '1';
fs = 100000;

%% get matlab generated results to compare

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

% generate index
index_corr = linspace(1,n_fft,n_fft);
index_corr_mat = linspace(1,n_fft-1,n_fft-1);

%% convert values to binary representation
% see page 45 of fft doc for values
if n_samples == 1024
    n_points_b = '000';
elseif n_samples == 2048
    n_points_b = '001';
elseif n_samples == 4096
    n_points_b = '010';
elseif n_samples == 8192
    n_points_b = '011';
elseif n_samples == 16384
    n_points_b = '100';
else
    n_points_b = '010';
end

if n_fft == 2048
    n_fft_b = '000';
elseif n_fft == 4096
    n_fft_b = '001';
elseif n_fft == 8192
    n_fft_b = '010';
elseif n_fft == 16384
    n_fft_b = '011';
else
    n_fft_b = '010';
end

cmd = bin2dec(strcat(run, ram_initialized, n_points_b, n_fft_b));

%% serial comms
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (8 * n_samples)); % 64 bits x n_fft points
set(s, 'OutputBufferSize', 1);
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');
fopen(s);
fwrite(s, cmd, 'uint8');
sdata = fread(s, (n_fft), 'int32');
cmd = bin2dec(strcat('0', ram_initialized, n_points_b, n_fft_b));
fwrite(s,cmd,'uint8');
fclose(s);

%% split real/imag parts, take magnitude
a = 1;
b = 1;
corr_fpga = zeros(n_fft,1);
sdata_im = zeros(1,n_fft);
sdata_re = zeros(1,n_fft);

for i = 1:2:(n_fft) 
    sdata_im(a) = sdata(i);
    a = a + 1;
end

for i = 2:2:(n_fft) 
    sdata_re(b) = sdata(i);
    b = b + 1;
end

% for i = 1:n_fft
%     corr_fpga(i,1) = int64(sqrt(sdata_im(1,i).^2 + sdata_re(1,i).^2));
% end   

% corr_fpga = sqrt(sdata_re.^2 + sdata_im.^2);
corr_fpga = (sdata_re + sdata_im);

%% plot quantized and non quantized xcorrs

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
% grid on;
% subplot(3,1,3);
% plot(index_corr, corr_fpga);
% grid on;
% title('FPGA xcorr');

%% plot result

figure;
subplot(3,1,1);
plot(index_corr_mat, corr_mat_quant);
title('Built in xcorr- quantized');
grid on;
subplot(3,1,2);
plot(index_corr, corr_man_quant);
title('Manual xcorr- quantized');
grid on;
subplot(3,1,3);
plot(index_corr, corr_fpga);
grid on;
title('FPGA xcorr');

