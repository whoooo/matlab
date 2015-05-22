%% set n_fft, run/rst, ram initialized, and sampling frequency (100kHz or 40kHz) 
n_samples = 2048;
n_fft = 4096;
run = '1'; % reset = 1 if run = 0
ram_initialized = '0';

%% convert values to binary representation
% see page 45 of fft doc for values
if n_samples == 1024
    n_points_b = '000';
elseif n_samples == 2048
    n_points_b = '001';
elseif n_samples == 4096
    n_points_b = '010';
elseif n_samples == 8192
    n_points_b = '011';
elseif n_samples == 16384
    n_points_b = '100';
else
    n_points_b = '010';
end

if n_fft == 2048
    n_fft_b = '000';
elseif n_fft == 4096
    n_fft_b = '001';
elseif n_fft == 8192
    n_fft_b = '010';
elseif n_fft == 16384
    n_fft_b = '011';
else
    n_fft_b = '010';
end

cmd = bin2dec(strcat(run, ram_initialized, n_points_b, n_fft_b));

%% serial comms
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'OutputBufferSize', 1);
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');
fopen(s);
fwrite(s, cmd, 'uint8');
% sdata = fread(s, (2 * n_points), 'int16');
fclose(s);