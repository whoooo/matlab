a=1;
b=1;
fs= 40000;
n_points = 2048;
fft_index = linspace(1, (fs), n_points);
n_index = linspace(1, n_points, n_points);
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (4 * n_points));
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');
fopen(s);
sdata = fread(s, (2 * n_points), 'int16');
fclose(s);
for i = 1:2:(2*n_points)
    if i == 1
        sdata_im(a) = sdata(i);
        a = a + 1;
    else
        sdata_im(a) = sdata(i);
        a = a + 1;
    end
end
for i = 2:2:(2*n_points)
    sdata_re(b) = sdata(i);
    b = b + 1;
end
for i = 1:n_points
    fft_mag(i,1) = int16(sqrt(sdata_im(1,i).^2 + sdata_re(1,i).^2));
end     
subplot(3,1,1);
plot(fft_index(1:length(fft_index)/2), sdata_re(1:length(sdata_re)/2));
title('FFT Real Values');
xlabel('Frequency');
ylabel('Magnitude');

subplot(3,1,2);
plot(fft_index(1:length(fft_index)/2), sdata_im(1:length(sdata_im)/2));
title('FFT Imaginary Values');
xlabel('Frequency');
ylabel('Magnitude');

subplot(3,1,3);
plot(fft_index(1:length(fft_index)/2), fft_mag(1:length(fft_mag)/2));
title('FFT Magnitude');
xlabel('Frequency');
ylabel('Magnitude');
