s3 = serial('COM6');
set(s2, 'BaudRate', 115200);
set(s2, 'InputBufferSize', 80000);
set(s2, 'Timeout', 10);
set(s2, 'ByteOrder', 'bigEndian');
fopen(s2);
xcdata = fread(s2, 40000, 'int16');
fclose(s2);