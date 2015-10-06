a = csvread('C:\Users\jtobin\Dropbox\sram_r_margin.csv', 1);

vin = a(:,1);
vout = a(:,2);

figure()
plot(vin, vout)
hold on 
grid on
plot(vout, vin)