% matlab automatically pads to next power of two size when doing corr

%% option 1: collect input data from adc
%n_samples = 40000;
% s2 = serial('COM6');
% set(s2, 'BaudRate', 115200);
% set(s2, 'InputBufferSize', (n_samples * 2));
% set(s2, 'Timeout', 15);
% set(s2, 'ByteOrder', 'bigEndian');
% fopen(s2);
% %xdata = fread(s2, n_samples, 'int16');
% xdata2 = fread(s2, n_samples, 'int16');
% fclose(s2);

%zeropad = transpose(linspace(0, 0, length(xdata))); %transpose not needed?
%xdatapad = vertcat(xdata, zeropad);
%xdata2pad = vertcat(xdata2, zeropad);
% indext = linspace(0, n_samples-1, n_samples);

% nfft = 2*2^nextpow2(n_samples);
% % indexc = linspace(1,(n_samples*2 - 1),(n_samples*2 - 1));
% % indext = linspace(0,1,n_samples);
% indexf = fs/2*linspace(0,1,nfft/2+1);
% indexcm = linspace(0, nfft, nfft);
% corrlength = length(xdata)+length(xdata2)-1;
% T = 1/fs;
% t = (0:n_samples-1)*T;

%% optione 2: generate input data
f1 = 5000;
f2 = 5000;
fs = 40000;
n_samples = 2048;
n_gen = linspace(1, n_samples, n_samples); % # of points to generate

%%% sine waves to cross correlate
xdata = sin(2*pi*f1/fs*n_gen);
xdata2 = sin(2*pi*f2/fs*n_gen);

%%% sine waves with half of the length zero padded (double original length)
zeropad = linspace(0, 0, n_samples);
xdata_p = cat(2, xdata, zeropad);
xdata2_p = cat(2, xdata2, zeropad);

%%% fft of each input
xdata_fft = abs(fft(xdata, n_samples));
xdata_p_fft = abs(fft(xdata_p, n_samples*2));
xdata2_fft = abs(fft(xdata2, n_samples));
xdata2_p_fft = abs(fft(xdata2_p, n_samples*2));

%%% unpadded indices for plots
indext = linspace(0, n_samples-1, n_samples);
indexf = fs/2*linspace(0, 1, n_samples/2);
%indexc_man = linspace(0, n_samples*2, (n_samples*2 - 1));
indexc_man = linspace(0, n_samples, (n_samples));
indexc_mat = linspace(0, n_samples*2, (n_samples*2 - 1));

%%% padded indices for plots
indext_p = linspace(0, (n_samples*2)-1, (2*n_samples));
indexf_p = fs/2*linspace(0, 1, n_samples);
% indexc_man_p = linspace(0, n_samples*4, (n_samples*4 - 1));
indexc_man_p = linspace(0, n_samples*2, (n_samples*2));
indexc_mat_p = linspace(0, n_samples*4, (n_samples*4 - 1));


%% use built in correlation function

%%% for unpadded data
corr_mat = xcorr(xdata, xdata2);

%%% for padded data
corr_mat_p = xcorr(xdata_p, xdata2_p);

%% manually correlate data

%%% for unpadded data
corr_man = fftshift(ifft(fft(xdata).*conj(fft(xdata2))));

%%% for padded data
corr_man_p = fftshift(ifft(fft(xdata_p).*conj(fft(xdata2_p))));

%% old- may still use
% xdata_fft = fft(xdata, nfft)/n_samples;
% xdata2_fft = fft(xdata2, nfft)/n_samples;
% xdata_fft = fft(xdatapad, nfft)/n_samples;
% xdata2_fft = fft(xdata2pad, nfft)/n_samples;
% corr_man = ifft(xdata_fft.*conj(xdata2_fft));

%% plot results

%%% plot xdata
subplot(5,2,1);
plot(indext, xdata);
grid minor;
legend('sample1');

%%% plot xdata2
subplot(5,2,2);
plot(indext, xdata2);
grid minor;
legend('sample2');

%%% plot f spectrum of xdata
subplot(5,2,3);
plot(indexf, xdata_fft(1:length(xdata_fft)/2));
grid minor;
legend('sample1 spectrum');

%%% plot f spectrum of xdata2
subplot(5,2,4);
plot(indexf, xdata2_fft(1:length(xdata_fft)/2));
grid minor;
legend('sample2 spectrum');

%%% plot f spectrum of padded xdata
subplot(5,2,5);
plot(indexf_p, xdata_p_fft(1:length(xdata_p_fft)/2));
grid minor;
legend('sample1(padded) spectrum');

%%% plot f spectrum of padded xdata2
subplot(5,2,6);
plot(indexf_p, xdata2_p_fft(1:length(xdata2_p_fft)/2));
grid minor;
legend('sample2(padded) spectrum');

%%% plot matlab unpadded correlation
subplot(5,2,7);
plot(indexc_mat, corr_mat);
grid minor;
legend('xcorr mat unpadded');

%%% plot matlab padded correlation
subplot(5,2,8);
plot(indexc_mat_p, corr_mat_p);
grid minor;
legend('xcorr mat padded');

%%% plot manual unpadded correlation
subplot(5,2,9);
plot(indexc_man, corr_man);
grid minor;
legend('xcorr man unpadded');

%%% plot manual padded correlation
subplot(5,2,10);
plot(indexc_man_p, corr_man_p);
grid minor;
legend('xcorr man padded');


