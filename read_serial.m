a=1;
b=1;
s = serial('COM6');
set(s, 'BaudRate', 115200);
set(s, 'InputBufferSize', 8192);
set(s, 'Timeout', 20);
set(s, 'ByteOrder', 'bigEndian');
fopen(s)
sdata = fread(s, [2048], 'int32');
fclose(s)
for i = 1:2:2048
    if i == 1
        sdata_im(a) = sdata(i);
        a = a + 1;
    else
        sdata_im(a) = sdata(i);
        a = a + 1;
    end
end
for i = 2:2:2048
    sdata_re(b) = sdata(i);
    b = b + 1;
end    
sdata_re_h = dec2hex(sdata_re, 8);
sdata_im_h = dec2hex(sdata_im, 8);
