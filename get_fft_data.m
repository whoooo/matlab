%% For use with adc_fft_to_pc2 Vivado project
% Must run once before taking samples

%% initialize parameters and get data
fs= 100000;
n_points = 4096;
ram_initialized = 1;
fft_index = linspace(1, (fs), n_points);
n_index = linspace(1, n_points, n_points);
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (2*2*n_points)); % 16 bit fft -> 2 bytes real + 2 bytes imag
set(s, 'OutputBufferSize', 1); % 1 byte command 
set(s, 'Timeout', 5);
set(s, 'ByteOrder', 'bigEndian');

%% convert values to binary representation
% see page 45 of fft doc for values
if n_points == 512
    n_points_b = 1001;
elseif n_points == 1024
    n_points_b = 1010;
elseif n_points == 2048
    n_points_b = 1011;
elseif n_points == 4096
    n_points_b = 1100;
elseif n_points == 8192
    n_points_b = 1101;
elseif n_points == 16384
    n_points_b = 1110;
elseif n_points == 32768
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


%% send command, then receive data
fopen(s);
fwrite(s, cmdrun, 'uint8');
sdata = fread(s, (2 * n_points), 'int16'); % 
fwrite(s, cmdrst, 'uint8');
fclose(s);

%% split real/imag parts, take magnitude
a = 1;
b = 1;
fft_mag = zeros(n_points,1);
fft_mag_shift = zeros(n_points,1);
sdata_im = zeros(1,n_points);
sdata_re = zeros(1,n_points);

for i = 1:2:(n_points*2) % *2
    sdata_im(a) = sdata(i);
    a = a + 1;
end

for i = 2:2:(n_points*2) % *2
    sdata_re(b) = sdata(i);
    b = b + 1;
end

fpga_fft = sdata_re + (sdata_im.*1j);
fpga_fft_conj = conj(fpga_fft);

for i = 1:n_points
    fft_mag(i,1) = sqrt(sdata_im(1,i).^2 + sdata_re(1,i).^2);
%     fft_mag_shift(i,1) = sqrt(sdata_im_shift(1,i).^2 + sdata_re_shift(1,i).^2);
end     

% sdata_re_shift = fftshift(sdata_re);
% sdata_im_shift = fftshift(sdata_im);

%% gunshot plots
figure
a= 2;
b= 1;

subplot(a,b,1);
plot(fft_index(1:length(fft_index)), (1/n_points) .* fft_mag(1:length(fft_mag)));
title('FFT Magnitude');
xlabel('Frequency');
ylabel('Magnitude');

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
