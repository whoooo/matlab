%% compare results of 16 bit input FFT and 32 bit input FFT using same input values (resize 16 bit to 32).
% use fft_16_32_compare Vivado project

% scaling_factor = 2^sum(scaling_sch);
% fft_mat = fft(input);
% mag_mat = abs(fft_mat)/scaling_factor;

% if zeropad = 1, max n_points to send is 4096
% if zeropad = 0, max n_points is be 8192

%% initialize values and set serial parameters 
n_points = 2048; % 2048 max
ram_init = 1; % use stored values from coe file if 1, otherwise use adc
fs = 40000;
zeropad = 0; % zero pad 2 * signal_length
run = 1;
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (n_points * (4*(1 + fft32)))); % take in 32 bits real/32 bits imaginary if fft32 = 1, else 16 bits each
set(s, 'OutputBufferSize', 1);
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');

%% convert command values to binary representation and send
% see page 45 of fft doc for values
if n_points == 512
    n_points_b = 1001;
elseif n_points == 1024
    n_points_b = 1010;
elseif n_points == 2048
    n_points_b = 1011;
elseif n_points == 4096
    n_points_b = 1100;
else
    n_points_b = 1010;
end

%                         7                   6             5          4               (3:0)                                
cmd = bin2dec(strcat(num2str(run), num2str(ram_init), num2str(zeropad), num2str(0), num2str(n_points_b)));
cmd_rst = bin2dec(strcat(num2str(0), num2str(ram_init), num2str(zeropad), num2str(0), num2str(n_points_b)));


%% write and read serial data

fopen(s);
fwrite(s, cmd, 'uint8');

if (fft32 == 1 && zero_pad == 1)
    data32pad = fread(s, (2 * n_points), 'int32');
elseif (fft32 == 1 && zero_pad == 0)
    data32 = fread(s, (2 * n_points), 'int32');
elseif (fft32 == 0 && zero_pad == 1)
    data16pad = fread(s, (2 * n_points), 'int16');
elseif (fft32 == 0 && zero_pad == 0)
    data16 = fread(s, (2 * n_points), 'int16');