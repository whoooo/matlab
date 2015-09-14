% Cross correlation based matching algorithm test using free firearm library sounds                        
% Uses fingerprints generated from gen_fingerprint

clear all

%% initialization

nfft = 4096; % actual number of samples taken will be 1/2 this, due to zero padding
fs = 48000; % desired sampling frequency
corrthresh = 8000;
add_noise = 0;

% gunshot files
fingerprint = load('Z:\jtobin\gunshots\fingerprintLib\f_domain\mat_files\R_27_s1_2048_4096_48k.txt');
fingerprint = fingerprint(:,1) + 1j.*fingerprint(:,2); % combine real and imaginary components 
sample = 'Z:\jtobin\gunshots\FreeFirearmLibrary\rawLibrary\I_15.wav'; %r27

% additional noise file
[z, fs_origZ] = audioread('Z:\jtobin\other_sounds\milan_traffic.wav');



%% Read and resample data. 

[y, fs_origY] = audioread(sample);

% (:,1) is first channel, (:,2) is second
% Downsample
ych1 = resample(y(:,1),fs,fs_origY);
ych2 = resample(y(:,2),fs,fs_origY);
ych1 = awgn(ych1, 10);

zch1 = resample(z(:,1), fs, fs_origZ);
zch2 = resample(z(:,2), fs, fs_origZ);

zch1 = zch1(1:length(ych1));
zch2 = zch2(1:length(ych1));


if add_noise == 1
    ych1 = normalize(ych1,1) + normalize(zch1, .6);
    ych2 = normalize(ych2,1) + normalize(zch2, .6);
end

% plot( linspace(1,length(zch1),length(zch1)), zch1);

%% threshold based correlation detection

corr_count1 = 1;
corr_count2 = 1;
zeropad = transpose(linspace(0, 0, nfft/2));

for k = 0 : 1 : (2*floor(length(ych1)/nfft) - 1) 
    
    sample1 = ych1((k*nfft/2 + 1):((k+1)*nfft/2));
    sample1pad = cat(1, sample1, zeropad);  
    sample1pad_f = fft(sample1pad, nfft);
    corr1 = ifft(sample1pad_f.*fingerprint); 
    if max(corr1) > corrthresh
        matches_index1(:,corr_count1) = (k*nfft/2 + 1) : ((k+1)*nfft/2); % time indices of matches
        matches_corr1(:,corr_count1) = corr1; % matching correlations
        index_corr1(:,corr_count1) = (k*nfft/2 + 1) : ((k+1)*nfft/2 + nfft/2); % indices for matching correlations
        matches_sample1(:,corr_count1) = sample1; % matching t data
        matches_fft1(:,corr_count1) = sample1pad_f; % matching f data
        corr_count1 = corr_count1 + 1;
    end;   
    
    sample2 = ych2((k*nfft/2 + 1):((k+1)*nfft/2));
    sample2pad = cat(1, sample2, zeropad);  
    sample2pad_f = fft(sample2pad,nfft);
    corr2 = ifft(sample2pad_f.*fingerprint); 
    if max(corr2) > corrthresh 
        matches_index2(:,corr_count2) = ((k*nfft/2 + 1):((k+1)*nfft/2));
        matches_corr2(:,corr_count2) =  corr2;
        index_corr2(:,corr_count2) = (k*nfft/2 + 1) : ((k+1)*nfft/2 + nfft/2); 
        matches_sample2(:,corr_count2) = sample2;
        matches_fft2(:,corr_count2) = sample2pad_f;
        corr_count2 = corr_count2 + 1;
    end; 

end;

%% generate indices
index_t = linspace(1,length(ych1),length(ych1));
index_f = transpose(fs/2*linspace(0,1,nfft));
index_fp = transpose(linspace(1,length(fingerprint),length(fingerprint)));


%% plot ch1 results
figure;
a = 3;
b = 2;

% ych1 time domain
subplot(a,b,[1,2]);
plot(index_t, ych1, 'b',...
    matches_index1, matches_sample1, 'r');
title('\bf Channel 1');
grid on;
xlim([1 length(ych1)])
xlabel('Index');
ylabel('Amplitude');

% correlation at match
subplot(a,b,[3,4]);
plot(index_corr1, matches_corr1);
title('\bf Correlation at match');
grid on;
xlim([1 length(ych1)])
xlabel('Correlation index');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,5);
plot(index_fp,ifft(conj(fingerprint), nfft));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% fingerprint f domain
subplot(a,b,6);
plot(index_f(1:floor(2*length(index_f)/5)), abs(normalize(fingerprint(1:floor(2*length(index_f)/5)),1)), 'b',...
    index_f(1:floor(2*length(index_f)/5)), abs(normalize(matches_fft1(1:floor(2*length(index_f)/5)),1)), 'r');
title('\bf Fingerprint/Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');

%% plot ch2 results
figure;
a = 3;
b = 2;

% ych2 time domain
subplot(a,b,[1,2]);
plot(index_t, ych2, 'b',...
    matches_index2, matches_sample2, 'r');
title('\bf Channel 2');
grid on;
xlim([1 length(ych2)])
xlabel('Index');
ylabel('Amplitude');

% correlation at match
subplot(a,b,[3,4]);
plot(index_corr2, matches_corr2);
title('\bf Correlation at match');
grid on;
xlim([1 length(ych2)])
xlabel('Correlation index');
ylabel('Amplitude');

% fingerprint t domain
subplot(a,b,5);
plot(index_fp,ifft(conj(fingerprint), nfft));
title('\bf Fingerprint');
grid on;
xlabel('Sample index');
ylabel('Amplitude');

% fingerprint f domain
subplot(a,b,6);
plot(index_f(1:floor(2*length(index_f)/5)), abs(normalize(fingerprint(1:floor(2*length(index_f)/5)),1)), 'b',...
    index_f(1:floor(2*length(index_f)/5)), abs(normalize(matches_fft2(1:floor(2*length(index_f)/5)),1)), 'r');
title('\bf Fingerprint/Match FFT');
grid on;
xlabel('Frequency');
ylabel('Amplitude');


%% plot ch2 results
% figure;
% a = 3;
% b = 2;
% 
% % ych2
% subplot(a,b,[1,2]);
% plot(index_t, ych1, 'b',...
%     matches_index2, matches_sample2, 'r');
% title('\bf Channel 2');
% grid on;
% xlim([1 length(ych2)])
% xlabel('Index');
% ylabel('Amplitude');
% 
% % fingerprint and match f domain
% subplot(a,b,[3,4]);
% plot(index_f(1:floor(2*length(index_f)/5)), abs(fingerprint(1:floor(2*length(index_f)/5))), 'b',...
%     index_f(1:floor(2*length(index_f)/5)), match_f2(1:floor(2*length(index_f)/5)), 'r');
% title('\bf Fingerprint/Match FFT');
% grid on;
% xlabel('Frequency');
% ylabel('Amplitude');
% 
% % fingerprint t domain
% subplot(a,b,5);
% plot(index_fp,ifft(conj(fingerprint), nfft));
% title('\bf Fingerprint');
% grid on;
% xlabel('Sample index');
% ylabel('Amplitude');
% 
% % correlation at match
% % subplot(a,b,6);
% % plot(index_corr2,match_corr2);
% % title('\bf Correlation at match');
% % grid on;
% % xlabel('Correlation index');
% % ylabel('Amplitude');


