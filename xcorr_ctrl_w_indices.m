clear all

%% set commands
n_samples = 2048;
n_fft = n_samples * 2;
run = '0';
run_once = '0';
rst = '0';
use_adc = '0';
fs = 48000;
send_results = '0';

%% convert command to binary representation

if (run == '1') && (rst == '1')
    run = '0';
    rst = '0';
    disp('Run and rst both set to 1');
end

% see page 45 of fft doc for values
if n_samples == 256
    n_points_b = '000';
elseif n_samples == 512
    n_points_b = '001';
elseif n_samples == 1024
    n_points_b = '010';
elseif n_samples == 2048
    n_points_b = '011';
elseif n_samples == 4096
    n_points_b = '100';
else
    n_points_b = '100';
end

cmd_config = bin2dec(strcat('0', '0', use_adc, run_once, n_points_b, send_results));
cmd = bin2dec(strcat(run, rst, use_adc, run_once, n_points_b, send_results));

%% serial comms
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', (8)); % 16 bit + 16 bit + 32 bit
set(s, 'OutputBufferSize', 1);
set(s, 'Timeout', 10);
set(s, 'ByteOrder', 'bigEndian');
fopen(s);
fwrite(s, cmd_config, 'uint8'); % send configuration command before sending run/rst command
fwrite(s, cmd, 'uint8');

index = 0;
n_detections = 0;
n_detections_total = 0;
event_match_indices = 0;
if run == '1' && send_results == '1' 
    index = index + 1;
    n_detections(index) = fread(s, 1, 'int16');
    n_detections_total(index) = fread(s, 1, 'int16'); 
    event_match_indices(index) = fread(s, 1, 'int32');       
end

fclose(s);