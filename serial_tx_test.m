%% set n_fft, run/rst, ram initialized, and sampling frequency (100kHz or 40kHz) 
n_points = 2048;
run = 1; % reset = 1 if run = 0
ram_initialized = 0;
fs = 40000; 

%% convert values to binary representation
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