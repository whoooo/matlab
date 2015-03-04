t = linspace(1, 8192, 8192);
f1 = 5000;
f2 = 5000;
fs = 40000;
sin5k = sin(2*pi*f1/fs*t);
sin7k = sin(2*pi*f2/fs*t);

quant_sin5k = int16(round(sin5k * 2^15));
quant_sin7k = int16(round(sin7k * 2^15));

%quant_sin1k_h = num2hex(quant_sin1k);
% for i = 1:1:1024
%     quant_sin1k_h() = sprintf('%x', typecast(int16(quant_sin1k(1, i)),'uint16'));
% end

subplot(2,1,1);
plot(t(1:512), quant_sin5k(1:512));
subplot(2,1,2);
plot(t(1:512), quant_sin7k(1:512));
