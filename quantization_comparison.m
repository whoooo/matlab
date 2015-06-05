% compare effects of normalization and quantization being performed before fft vs after fft
n_samples = 2048;
n_fft = 4096;
fs = 100000;

%% initialization

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


%% quantize signal, take fft, compare quantized version of ifft to non quantized

y1_quant = double(int16(y1shock/max(abs(y1shock)).*32768));
y1_quant_fft = fft(y1_quant);
y1_quant_fft_quant = double(int16(y1_quant_fft/max(abs(y1_quant_fft)).*32768));
y1_quant_ifft = int16(ifft(y1_quant_fft_quant));


plot_section1 =1;

if plot_section1 == 1
    
figure;
subplot(3,1,1);
plot(linspace(0,length(y1shock),length(y1shock)), y1shock);
title('Y1shock, no quant');
subplot(3,1,2);
plot(linspace(0,length(y1_quant),length(y1_quant)), y1_quant);
title('Y1shock, quant');
subplot(3,1,3);
plot(linspace(0,length(y1_quant),length(y1_quant)), y1_quant_ifft);
title('IFFT of quantized FFT data');

figure;
subplot(2,1,1);
plot(linspace(0,length(y1_quant_fft),length(y1_quant_fft)), abs(y1_quant_fft));
title('fft');
subplot(2,1,2);
plot(linspace(0,length(y1_quant_fft),length(y1_quant_fft)), abs(y1_quant_fft_quant));
title('quantized fft');

end;




%% quantize both time signals before taking fft, do not quantize fft results

% quantize and normalize time domain signal
% must convert to double after int16 as fft will only work on uint16, not
y1s_quant = double(int16(y1shock/abs(max(y1shock)).*32768));
y2s_quant = double(int16(y2shock/abs(max(y2shock)).*32768));

% zero pad sample to at least twice its length
zeropad = transpose(linspace(0, 0, n_fft - n_samples));
y1shockpad_quant = cat(1,y1s_quant,zeropad);
y2shockpad_quant = cat(1,y2s_quant,zeropad);

% take manual cross correlation
corr_man_y1t = ifft(fft(y1shockpad_quant, n_fft).*conj(fft(y2shockpad_quant, n_fft)));
corr_man_y2t = ifft(conj(fft(y1shockpad_quant, n_fft)).*(fft(y2shockpad_quant, n_fft)));

% use built in xcorr function
corr_mat = xcorr(y1s_quant, y2s_quant);

% generate index
index_corr = linspace(1,n_fft,n_fft);
index_corr_mat = linspace(1,n_fft-1,n_fft-1);

plot_section2 = 0;

if plot_section2 == 1
    
figure;
subplot(3,1,1)
plot(index_corr, fftshift(corr_man_y1t));
title('y1t manual');
subplot(3,1,2)
plot(index_corr, fftshift(corr_man_y2t));
title('y2t manual');
subplot(3,1,3)
plot(index_corr_mat, corr_mat)
title('Built in xcorr');

end;

%% quantize one time signal before taking fft, quantize results of fft

% % take complex conjugate of fft
% y1s_f_quant = conj(fft(y1shockpad_quant, nfft));
% y2s_f_quant = conj(fft(y2shockpad_quant, nfft));
% 
% % normalize and quantize the complex conjugate of fft data
% y1s_fp_norm = y1s_fp_f/abs(max(y1s_fp_f));
% y1s_fp_quant = int16(y1s_fp_norm .*32768);
% y2s_fp_norm = y2s_fp_f/abs(max(y2s_fp_f));
% y2s_fp_quant = int16(y2s_fp_norm .*32768);