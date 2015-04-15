% matlab automatically pads to next power of two size when doing corr
% scaling_factor = 2^sum(scaling_sch);
% fft_mat = fft(input);
% mag_mat = abs(fft_mat)/scaling_factor;

%% initialize values and set serial parameters 
n_samples = 8192;
n_points = n_samples;
fs = 100000;
run = 1;
ram_initialized = 0;
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (n_samples * 2));
set(s, 'OutputBufferSize', 1);
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');

%% generate indices
indext = linspace(0, n_samples-1, n_samples); % time index
indexf = fs/2*linspace(0, 1, n_samples/2); % frequency index
indexf_p = fs/2*linspace(0, 1, n_samples); % padded frequency index
indexc_man_p = linspace(0, (n_samples*2 - 1), (n_samples*2)); % manual correlation index
indexc_mat = linspace(0, (n_samples*2 - 1), (n_samples*2 - 1)); % built in correlation function index

%% convert command values to binary representation
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

cmd = bin2dec(strcat(num2str(run), num2str(fs_b), num2str(ram_initialized), num2str(0), num2str(n_points_b)));
cmd_rst = bin2dec(strcat(num2str(0), num2str(fs_b), num2str(ram_initialized), num2str(0), num2str(n_points_b)));

%% option 1: collect input data from adc
% fopen(s);
% fwrite(s, cmd, 'uint8');
% % xdata = fread(s, n_samples, 'int16');
% xdata2 = fread(s, n_samples, 'int16');
% fwrite(s, cmd_rst, 'uint8');
% fclose(s);
% % xdata =xdata2;

%% option 2: generate input data
f1 = 1000;
f2 = 1000;
n_gen = linspace(1, n_samples, n_samples); % # of points to generate

%%% generate sine waves with frequencies f1 and f2
xdata = transpose(32768*sin(2*pi*f1/fs*n_gen));
xdata2 = transpose(32768*sin(2*pi*f2/fs*n_gen));

% %%% generate rectangular function
% xdata = linspace(1, 1, n_samples);
% xdata2 = linspace(1, 1, n_samples);

%% zero pad signals
zeropad = transpose(linspace(0, 0, n_samples));
xdatapad = cat(1, xdata, zeropad);
xdata2pad = cat(1, xdata2, zeropad);

%% scale input samples
% xdata = xdata * 10/32768;
% xdata2 = xdata2 * 10/32768;

%% fft of signals
xdata_fft = abs(fft(xdata, n_samples));
xdatapad_fft = abs(fft(xdatapad, n_samples*2));
xdata2_fft = abs(fft(xdata2, n_samples));
xdata2pad_fft = abs(fft(xdata2pad, n_samples*2));

%% correlate data

%%% using built in matlab function
corr_mat = xcorr(xdata, xdata2);

%%% manual correlation using padded data
corr_man_p = fftshift(ifft(fft(xdatapad).*conj(fft(xdata2pad))));

%% plot results
figure
%%% plot xdata
subplot(4,2,1);
plot(indext, xdata);
grid minor;
legend('sample1');
xlabel('Sample Index'), ylabel('Magnitude')
axis([0 n_samples -32768 32768])

%%% plot xdata2
subplot(4,2,2);
plot(indext, xdata2);
grid minor;
legend('sample2');
xlabel('Sample Index'), ylabel('Magnitude')
axis([0 n_samples -32768 32768])

%%% plot f spectrum of xdata
subplot(4,2,3);
plot(indexf, xdata_fft(1:length(xdata_fft)/2));
grid minor;
xlabel('Frequency'), ylabel('Magnitude')
legend('sample1 spectrum');

%%% plot f spectrum of xdata2
subplot(4,2,4);
plot(indexf, xdata2_fft(1:length(xdata_fft)/2));
grid minor;
xlabel('Frequency'), ylabel('Magnitude')
legend('sample2 spectrum');

%%% plot f spectrum of padded xdata
subplot(4,2,5);
plot(indexf_p, xdatapad_fft(1:length(xdatapad_fft)/2));
grid minor;
xlabel('Frequency'), ylabel('Magnitude')
legend('sample1 padded spectrum');

%%% plot f spectrum of padded xdata2
subplot(4,2,6);
plot(indexf_p, xdata2pad_fft(1:length(xdata2pad_fft)/2));
grid minor;
xlabel('Frequency'), ylabel('Magnitude')
legend('sample2 padded spectrum');

%%% plot matlab unpadded correlation
subplot(4,2,7);
plot(indexc_mat, corr_mat);
grid minor;
xlabel('Index'), ylabel('Magnitude')
legend('xcorr built in');

%%% plot manual padded correlation
subplot(4,2,8);
plot(indexc_man_p, corr_man_p);
grid minor;
xlabel('Index'), ylabel('Magnitude')
legend('xcorr manual');


