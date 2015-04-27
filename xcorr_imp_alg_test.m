% Cross correlation based matching algorithm test using free firearm library sounds                        
% Uses fingerprints generated in gunshot_audio_read.m

% clear all

%% initialization

nfft = 4096; % actual number of samples taken will be 1/2 this, due to zero padding
fs = 100000; % desired sampling frequency
fingerprint = load('Z:\jtobin\gunshots\fingerprintLib\R_27_s1_1024_4096_100k.mat');
fingerprint_half = load('Z:\jtobin\gunshots\fingerprintLib\R_27_s1_1024_4096_100k_half_freq.mat');
sample = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\R_27.wav'; %s11

fingerprint = fingerprint.y1s_fp; 
fingerprint_half = fingerprint_half.y1s_fp;


%% Read and resample data. 

[y, fs_orig] = audioread(sample);

% (:,1) is first channel, (:,2) is second
ych1 = y(:,1);
ych2 = y(:,2);

% Downsample
ych1 = resample(ych1,fs,fs_orig);
ych2 = resample(ych2,fs,fs_orig);

%% "slide" sample across window and check correlation
% correlation = fftshift(ifft(fft(xdatapad).*conj(fft(xdata2pad))));
% add threshold and count number of detections?
% sample length is half of nfft length to allow for zero padding
max1 = 0;
max2 = 0;
ch1count = 0;
ch2count = 0;
corrthresh = 5;
zeropad = transpose(linspace(0, 0, nfft/2));

for i = 0 : 1 : (floor(length(ych1)/nfft) - 1) 
    
    sample1 = ych1((i*nfft/2 + 1):((i+1)*nfft/2));
    sample1pad = cat(1, sample1, zeropad);  
%     corr1 = xcorr(sample1, y1shock);
    corr1 = ifft(fft(sample1pad,nfft).*fingerprint);   
    if max(corr1) > max1 
        max1 = max(corr1);
        match_sample1 = sample1;
        match_corr1 = corr1;
        match_i1 = i;
    end;   
    
    sample2 = ych2((i*nfft/2 + 1):((i+1)*nfft/2));
    sample2pad = cat(1, sample2, zeropad);  
%     corr2 = xcorr(sample2, y2shock);
    corr2 = ifft(fft(sample2pad,nfft).*fingerprint); 
    if max(corr2) > max2 
        max2 = max(corr2);
        match_sample2 = sample2;
        match_corr2 = corr2;
        match_i2 = i;
    end; 

end;

%% "slide" sample across window and check correlation using half spectrum
% correlation = fftshift(ifft(fft(xdatapad).*conj(fft(xdata2pad))));
% add threshold and count number of detections?
% sample length is half of nfft length to allow for zero padding
max_half = 0;
zeropad_half = transpose(linspace(0, 0, nfft/2));

for i = 0 : 1 : (floor(length(ych1)/nfft) - 1) 
    
    sample1 = ych1((i*nfft/2 + 1):((i+1)*nfft/2));
    sample1pad = cat(1, sample1, zeropad);
    samp1_f = fft(sample1pad,nfft);
    samp1_f_half = samp1_f(1:2048);
    corr1_half = ifft(samp1_f_half.*fingerprint_half);   
    if max(corr1_half) > max_half
        max_half = max(corr1_half);
        match_sample1_half = sample1;
        match_corr1_half = corr1_half;
        match_i1_half = i;
    end;   
    
end;

%% take fft of sample to inspect f domain
match_f1 = abs(fft(match_sample1, nfft));
match_f2 = abs(fft(match_sample2, nfft));
match_half = abs(fft(match_sample1_half, nfft));

%% generate indices
index_t1 = transpose((1000/fs).*linspace(1,length(ych1),length(ych1)));
index_t2 = transpose((1000/fs).*linspace(1,length(ych1),length(ych1)));
index_fp_f = transpose(fs/2*linspace(0,1,nfft));
index_match1 = transpose((1000/fs).*linspace(match_i1*nfft/2 + 1, (match_i1 + 1)*nfft/2, nfft/2));
index_match2 = transpose((1000/fs).*linspace(match_i2*nfft/2 + 1, (match_i2 + 1)*nfft/2, nfft/2));
index_match_half = transpose((1000/fs).*linspace(match_i1_half*nfft/2 + 1, (match_i1_half + 1)*nfft/2, nfft/2));
index_corr1 = transpose(linspace(1,length(match_corr1),length(match_corr1)));
index_corr2 = transpose(linspace(1,length(match_corr1),length(match_corr1)));
index_corr_half = transpose(linspace(1,length(match_corr1_half),length(match_corr1_half)));
index_fp = transpose(linspace(1,length(fingerprint),length(fingerprint)));


%% plot ch1 results
figure;
a = 3;
b = 2;

% ych1
subplot(a,b,[1,2]);
plot(index_t1,ych1, 'b', index_match1, match_sample1, 'r');
title('\bf Channel 1');
grid on;
xlabel('Index');
ylabel('Amplitude');

% fingerprint f domain
subplot(a,b,[3,4]);
% plot(index_fp_f,abs(fingerprint(1:length(fingerprint)/2)));
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), abs(fingerprint(1:floor(2*length(index_fp_f)/5))), 'b', index_fp_f(1:floor(2*length(index_fp_f)/5)), match_f1(1:floor(2*length(index_fp_f)/5)), 'r');
% plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), abs(fingerprint(1:floor(2*length(index_fp_f)/5))));
title('\bf Fingerprint/Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,5);
plot(index_fp,ifft(conj(fingerprint), nfft));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% correlation at match
subplot(a,b,6);
plot(index_corr1,match_corr1);
title('\bf Correlation at match');
grid on;
xlabel('Correlation index');
ylabel('Amplitude');



%% plot ch2 results
figure;
a = 3;
b = 2;

% ych2
subplot(a,b,[1,2]);
plot(index_t1,ych2, 'b', index_match2, match_sample2, 'r');
title('\bf Channel 2');
grid on;
xlabel('Index');
ylabel('Amplitude');

% fingerprint and match f domain
subplot(a,b,[3,4]);
% plot(index_fp_f,abs(fingerprint(1:length(fingerprint)/2)));
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), abs(fingerprint(1:floor(2*length(index_fp_f)/5))), 'b', index_fp_f(1:floor(2*length(index_fp_f)/5)), match_f2(1:floor(2*length(index_fp_f)/5)), 'r');
title('\bf Fingerprint/Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,5);
plot(index_fp,ifft(conj(fingerprint), nfft));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% correlation at match
subplot(a,b,6);
plot(index_corr2,match_corr2);
title('\bf Correlation at match');
grid on;
xlabel('Correlation index');
ylabel('Amplitude');

%% plot half spectrum results

figure;
a = 3;
b = 2;

% ych1
subplot(a,b,[1,2]);
plot(index_t1,ych1, 'b', index_match_half, match_sample1_half, 'r');
title('\bf Channel 1');
grid on;
xlabel('Index');
ylabel('Amplitude');

% fingerprint and match f domain
subplot(a,b,[3,4]);
% plot(index_fp_f,abs(fingerprint(1:length(fingerprint)/2)));
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), abs(fingerprint(1:floor(2*length(index_fp_f)/5))), 'b', index_fp_f(1:floor(2*length(index_fp_f)/5)), match_half(1:floor(2*length(index_fp_f)/5)), 'r');
title('\bf Half Spectrum Fingerprint/Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,5);
plot(index_fp,ifft(conj(fingerprint), nfft));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% correlation at match
subplot(a,b,6);
plot(index_corr_half,match_corr1_half);
title('\bf Half Spectrum Correlation at match');
grid on;
xlabel('Correlation index');
ylabel('Amplitude');

