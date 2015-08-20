%% For use with adc_fft_to_pc2 Vivado project
% Must run once before taking samples

clear all

%% initialize parameters and get data
fs= 100000;
n_fft = 4096;
n_samples = 2048;
ram_initialized = 1;
fft_index = linspace(1, (fs), n_fft);
n_index = linspace(1, n_fft, n_fft);
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (2*2*n_fft)); % 16 bit fft -> 2 bytes real + 2 bytes imag
set(s, 'OutputBufferSize', 1); % 1 byte command 
set(s, 'Timeout', 5);
set(s, 'ByteOrder', 'bigEndian');

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

% fft of original signal
y1s_fft = fft(y1shockpad, n_fft);
y2s_fft = fft(y2shockpad, n_fft);

% fft of quantized signal
y1s_quant_fft = fft(y1shockpad_quant, n_fft);
y2s_quant_fft = fft(y2shockpad_quant, n_fft);

%% convert command values to binary representation
% see page 45 of fft doc for values
if n_fft == 512
    n_points_b = 1001;
elseif n_fft == 1024
    n_points_b = 1010;
elseif n_fft == 2048
    n_points_b = 1011;
elseif n_fft == 4096
    n_points_b = 1100;
elseif n_fft == 8192
    n_points_b = 1101;
elseif n_fft == 16384
    n_points_b = 1110;
elseif n_fft == 32768
    n_points_b = 1111;
else
    n_points_b = 1010;
end

% fs = 40kHz if 0, 100kHz if 1
if fs == 40000
    fs_b = 0;
else
    fs_b = 1;   
end

cmdrun = bin2dec(strcat(num2str(1), num2str(fs_b), num2str(ram_initialized), num2str(0), num2str(n_points_b)));
cmdrst = bin2dec(strcat(num2str(0), num2str(fs_b), num2str(ram_initialized), num2str(0), num2str(n_points_b)));


%% send command and wait for data
fopen(s);
fwrite(s, cmdrun, 'uint8');
sdata = fread(s, (2 * n_fft), 'int16'); % 
fwrite(s, cmdrst, 'uint8');
fclose(s);

%% split real/imag parts, take magnitude
a = 1;
b = 1;
fft_mag = zeros(n_fft,1);
fft_mag_shift = zeros(n_fft,1);
sdata_im = zeros(1,n_fft);
sdata_re = zeros(1,n_fft);

for i = 1:2:(n_fft*2) % *2
    sdata_im(a) = sdata(i);
    a = a + 1;
end

for i = 2:2:(n_fft*2) % *2
    sdata_re(b) = sdata(i);
    b = b + 1;
end

fpga_fft = sdata_re + (sdata_im.*1j);
fpga_fft_conj = conj(fpga_fft);

% for i = 1:n_points
%     fft_mag(i,1) = sqrt(sdata_im(1,i).^2 + sdata_re(1,i).^2);
% %     fft_mag_shift(i,1) = sqrt(sdata_im_shift(1,i).^2 + sdata_re_shift(1,i).^2);
% end     

% fft_mag = abs(fpga_fft);

% sdata_re_shift = fftshift(sdata_re);
% sdata_im_shift = fftshift(sdata_im);

%% gunshot plots
figure
a= 3;
b= 1;

% subplot(a,b,1);
% plot(fft_index(1:length(fft_index)/2), (1/n_points) .* abs(fpga_fft(1:length(fpga_fft)/2)), 'g', ...
%     fft_index(1:length(fft_index)/2), abs(y1s_quant_fft(1:length(y1s_quant_fft)/2)), 'b' );
% title('FFT Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');

scale = 8192;

subplot(a,b,1);
plot(fft_index(1:length(fft_mag)/5), abs(fpga_fft(1:length(fft_mag)/5)))
title('FPGA FFT');
subplot(a,b,2);
plot(fft_index(1:length(fft_index)/5), abs(y1s_quant_fft(1:length(y1s_quant_fft)/5))/scale );
title('Matlab Quantized FFT');
subplot(a,b,3);
plot(fft_index(1:length(fft_index)/5), abs(y1s_fft(1:length(y1s_fft)/5))/scale );
title('Matlab FFT');

figure 

plot(fft_index(1:length(fft_mag)/5), abs(fpga_fft(1:length(fft_mag)/5)), 'r' , fft_index(1:length(fft_index)/5), abs(y1s_quant_fft(1:length(y1s_quant_fft)/5))/scale*2 , 'b');
legend('FPGA FFT', 'Matlab FFT')

% figure 
% plot(fft_index(1:length(fft_index)/2), abs(y1s_quant_fft(1:length(y1s_quant_fft)/2)));

% subplot(a,b,2);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* fft_mag_shift(1:length(fft_mag_shift)));
% title('Shifted FFT Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');


%% plots
% a = 3;
% b = 1;
% 
% subplot(a,b,1);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* sdata_re(1:length(sdata_re)));
% title('FFT Real Values');
% xlabel('Frequency');
% ylabel('Magnitude');

% subplot(a,b,2);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* sdata_re_shift(1:length(sdata_re_shift)));
% title('Shifted FFT Real Values');
% xlabel('Frequency');
% ylabel('Magnitude');

% subplot(a,b,2);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* sdata_im(1:length(sdata_im)));
% title('FFT Imaginary Values');
% xlabel('Frequency');
% ylabel('Magnitude');

% subplot(a,b,4);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* sdata_im_shift(1:length(sdata_im_shift)));
% title('Shifted FFT Imaginary Values');
% xlabel('Frequency');
% ylabel('Magnitude');

% subplot(a,b,3);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* fft_mag(1:length(fft_mag)));
% title('FFT Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');

% subplot(a,b,6);
% plot(fft_index(1:length(fft_index)), (1/n_points) .* fft_mag_shift(1:length(fft_mag_shift)));
% title('Shifted FFT Magnitude');
% xlabel('Frequency');
% ylabel('Magnitude');
