nfft = 4096; 
fs = 100000;
sample = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\R_27.wav'; 

% read audio and save channel 1
[y, fs_orig] = audioread(sample);
ych1 = y(:,1);

% Downsample
ych1 = resample(ych1,fs,fs_orig);
y1shock = ych1(44000:(44000+nfft-1));

% create windows
wb = window(@blackman, nfft);
wh = window(@hamming, nfft);
wg = window(@gausswin, nfft);

% apply windows
ych1_wb = y1shock.*wb;
ych1_wh = y1shock.*wh;
ych1_wg = y1shock.*wg;

% take FFTs
fft_rect = abs(fft(y1shock, nfft));
fft_wb = abs(fft(ych1_wb, nfft));
fft_wh = abs(fft(ych1_wh, nfft));
fft_wg = abs(fft(ych1_wg, nfft));

% create indices
index_f = transpose(fs/2*linspace(0,1,nfft));
index_t = transpose(linspace(1, nfft, nfft));

% plot
figure
subplot(2,1,1)
plot(index_t, y1shock);
xlabel('Index')
ylabel('Amplitude')
subplot(2,1,2)
plot(index_f(1:length(index_f)/2), fft_orig(1:length(fft_orig)/2), 'b'...
,index_f(1:length(index_f)/2), fft_wb(1:length(fft_orig)/2), 'r' ...
,index_f(1:length(index_f)/2), fft_wh(1:length(fft_orig)/2), 'm' ...
,index_f(1:length(index_f)/2), fft_wg(1:length(fft_orig)/2), 'c')
xlabel('Frequency')
ylabel('Amplitude')
legend('Rectangular', 'Blackman', 'Hamming', 'Gaussian')
grid on



