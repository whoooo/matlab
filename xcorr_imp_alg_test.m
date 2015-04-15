% Cross correlation based matching algorithm test using free firearm library sounds                        
% Uses fingerprints generated in gunshot_audio_read.m

clear all

%% initialization

nfft = 4096;
corrthresh = 1;
fs = 100000; % desired sampling frequency
fingerprint = load('Z:\jtobin\gunshots\fingerprintLib\R_27_s2_1024_4096_100k.mat');
sample = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\R_27.wav'; %s11

fingerprint = fingerprint.y2s_fp; 


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
max1 = 0;
max2 = 0;
zeropad = transpose(linspace(0, 0, nfft));

for i = 0 : 0.5 : (floor(length(ych1)/nfft) - 1) % change if using different samples instead of different channels
    sample1 = ych1((i*nfft + 1):((i+1)*nfft));
    sample1pad = cat(1, sample1, zeropad);  
    corr1 = fftshift(ifft(fft(sample1pad)).*fingerprint);   
    if max(corr1) > max1 
        max1 = max(corr1);
        match_sample1 = sample1;
        match_corr1 = corr1;
        match_i1 = i;
    end;   
    
    sample2 = ych2( (i*nfft + 1):((i+1)*nfft) );
    sample2pad = cat(1, sample2, zeropad);  
    corr2 = fftshift(ifft(fft(sample2pad, 2*nfft)).*fingerprint); 
    if max(corr2) > max2 
        max2 = max(corr2);
        match_sample2 = sample2;
        match_corr2 = corr2;
        match_i2 = i;
    end; 

end;

%% take fft of sample to inspect f domain
match_f1 = abs(fft(match_sample1, nfft));
match_f2 = abs(fft(match_sample2, nfft));

%% generate indices
index_t1 = transpose((1000/fs).*linspace(1,length(ych1),length(ych1)));
index_t2 = transpose((1000/fs).*linspace(1,length(ych1),length(ych1)));
index_fp_f = transpose(fs/2*linspace(0,1,nfft));
index_match1 = transpose((1000/fs).*linspace(match_i1*nfft + 1, (match_i1 + 1)*nfft, nfft));
index_match2 = transpose((1000/fs).*linspace(match_i2*nfft + 1, (match_i2 + 1)*nfft, nfft));
index_corr1 = transpose(linspace(1,length(match_corr1),length(match_corr1)));
index_corr2 = transpose(linspace(1,length(match_corr1),length(match_corr1)));
index_fp = transpose(linspace(1,length(fingerprint),length(fingerprint)));


%% plot ch1 results
figure;
a = 3;
b = 2;

% ych1
subplot(a,b,1);
plot(index_t1,ych1);
title('\bf Channel 1');
grid on;
xlabel('Time (mS)');
ylabel('Amplitude');

% correlation at match
subplot(a,b,2);
plot(index_corr1,match_corr1);
title('\bf Correlation at match');
grid on;
xlabel('Correlation index');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,3);
plot(index_fp,ifft(conj(fingerprint), nfft*2));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% fingerprint f domain
subplot(a,b,4);
% plot(index_fp_f,abs(fingerprint(1:length(fingerprint)/2)));
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), abs(fingerprint(1:floor(2*length(index_fp_f)/5))));
title('\bf Fingerprint FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');

% match
subplot(a,b,5);
plot(index_match1,match_sample1);
title('\bf Match');
grid on;
xlabel('Time (mS)');
ylabel('Amplitude');

% match fft
subplot(a,b,6);
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), match_f1(1:floor(2*length(index_fp_f)/5)));
title('\bf Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');


%% plot ch2 results
figure;
a = 3;
b = 2;

% ych2
subplot(a,b,1);
plot(index_t2,ych2);
title('\bf Channel 2');
grid on;
xlabel('Time (mS)');
ylabel('Amplitude');

% correlation at match
subplot(a,b,2);
plot(index_corr2,match_corr2);
title('\bf Correlation at match');
grid on;
xlabel('Correlation index');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,3);
plot(index_fp,ifft(conj(fingerprint), nfft*2));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% fingerprint f domain
subplot(a,b,4);
% plot(index_fp_f,abs(fingerprint(1:length(fingerprint)/2)));
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), abs(fingerprint(1:floor(2*length(index_fp_f)/5))));
title('\bf Fingerprint FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');

% match
subplot(a,b,5);
plot(index_match2,match_sample2);
title('\bf Match');
grid on;
xlabel('Time (mS)');
ylabel('Amplitude');

% match fft
subplot(a,b,6);
plot(index_fp_f(1:floor(2*length(index_fp_f)/5)), match_f2(1:floor(2*length(index_fp_f)/5)));
title('\bf Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');


