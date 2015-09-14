clear all
% 1024 | 2048 |4096 | 8192
samp_bram = [0.5 0.5 1 2];
sampf_bram = [1 2 4 7.5];
fp_bram = [1 2 4 7.5];
fft_bram = [7 7 11 22];
ifft_bram = [14 14 22 44];
n_fp = 1:120;
line = linspace(135,135,length(n_fp));

total_bram = [ 0 0 0 0];

for i = 1:4
    for j = 1:length(n_fp)
    total_bram(j,i) = samp_bram(i) +fft_bram(i) + sampf_bram(i) + ifft_bram(i) + fp_bram(i)*n_fp(j);
    end
end

total_bram(:,5) = line;

plot(n_fp, total_bram)
xlabel('Number of Fingerprints');
ylabel('BRAM Used')
legend('NFFT = 1024','NFFT = 2048','NFFT = 4096','NFFT = 8192')
xlim([ 1 120])
grid on