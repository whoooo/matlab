%% generate quantized sine waves of specific length, frequency, and sample frequency

sinlength = 2048;
t = linspace(1, sinlength, sinlength);
tlong = linspace(1,sinlength*2,sinlength*2);
% f_incr = linspace(1, 10, sinlength);
f1 = 5000;
f2 = 5000;
fs = 40000;
sin1 = sin(2*pi*f1/fs*t);
sin2 = sin(2*pi*f2/fs*t);
% sin_sum = cat(2,sin1,sin2);

quant_sin1 = int8(round(sin1 * 2^7));
quant_sin2 = int16(round(sin2 * 2^15));
% quant_sin_sum = int16(round(sin_sum * 2^15));

% xdata_fft = fft(sin_sum, 2048);
% indexf = fs/2*linspace(0, 1, sinlength*2);

subplot(2,1,1);
plot(t(1:512), quant_sin1(1:512));
subplot(2,1,2);
plot(t(1:512), quant_sin2(1:512));

% plot(tlong, sin_sum);
% plot(indexf, abs(xdata_fft(1:length(xdata_fft)/2)));
