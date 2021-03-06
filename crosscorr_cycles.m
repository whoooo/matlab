clear all

nfft = [1024 2048 4096 8191];
n_fingerprints = [1:250];
% n_fingerprints = 10;
fftcycles = [3458 7321 14489 30896];
clockcycles = [0 0 0 0];

for i = 1 : 4
    for j = 1 : 250
     clockcycles(j,i) = 1*nfft(i) + fftcycles(i) + n_fingerprints(j)*( 1*nfft(i) + fftcycles(i) );
    end
end

plot(n_fingerprints, clockcycles./100000)
xlabel('Number of Fingerprints');
ylabel('Total Run Time (mS)')
legend('NFFT = 1024','NFFT = 2048','NFFT = 4096','NFFT = 8192')
xlim([ 1 250])
grid on
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)


sampling_f = [ 48 100 ];
sample_period =[ 0 0];
for i = 1:4
    for j = 1:2
    sample_period(j,i) = nfft(i)/sampling_f(j);
    end
end



