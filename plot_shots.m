clear all

fdes = 48000; % desired sampling frequency
n_fft = 4096;
samplelength = 70000;

%% open files
x_39 = transpose(load('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\X_39_s1_2048_48k.txt'));
q_30 = transpose(load('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\Q_30_s1_2048_48k.txt'));
f_46 = transpose(load('Z:\jtobin\gunshots\fingerprintLib\time_domain\mat_files\F_46_s1_2048_48k.txt'));

lc9 = 'Z:\jtobin\gunshots\lc9.wav';
[lc9, fs] = audioread(lc9);
% lc9 = resample(lc9,fdes,fs);
lcp = 'Z:\jtobin\gunshots\lcp.wav';
[lcp, fs] = audioread(lcp); 
% lcp = resample(lcp,fdes,fs);
twenty2 = 'Z:\jtobin\gunshots\twenty2.wav';
[twenty2, fs] = audioread(twenty2);
% twenty2 = resample(twenty2,fdes,fs);

%% create indices
index_t = transpose(linspace(1,length(x_39),length(x_39)));
index_f = transpose(48000/2*linspace(0,1,n_fft/2));
index_t2 = transpose(linspace(1,length(lc9),length(lc9)));
index_t3 = transpose(linspace(1,length(lcp),length(lcp)));
index_t4 = transpose(linspace(1,length(twenty2),length(twenty2)));

%% cut recordings down
shockstartlc9 = 20000; 
lc9shock = lc9(shockstartlc9:(shockstartlc9+samplelength-1));
shockstartlcp = 10000;
lcpshock = lcp(shockstartlcp:(shockstartlcp+samplelength-1));
shockstarttwenty2 = 20000;
twenty2shock = twenty2(shockstarttwenty2:(shockstarttwenty2+samplelength-1));

index_shock = linspace(0,length(lc9shock),length(lc9shock));

%% plot recordings

subplot(3,1,1)
plot(index_t2, lc9, 'b', index_t2(shockstartlc9:(shockstartlc9+samplelength)), lc9(shockstartlc9:(shockstartlc9+samplelength)), 'r');
title('LC9');
subplot(3,1,2)
plot(index_t3, lcp, 'b', index_t3(shockstartlcp:(shockstartlcp+samplelength)), lcp(shockstartlcp:(shockstartlcp+samplelength)), 'r');
title('LCP');
subplot(3,1,3)
plot(index_t4, twenty2, 'b', index_t4(shockstarttwenty2:(shockstarttwenty2+samplelength)), twenty2(shockstarttwenty2:(shockstarttwenty2+samplelength)), 'r');
title('.22');

set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

%% plot cut recordings
figure
plot(index_shock, lc9shock);
xlabel('Time');
ylabel('Amplitude');
title('LC9')
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

figure
plot(index_shock, lcpshock);
xlabel('Time');
ylabel('Amplitude');
title('LCP')
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

figure
plot(index_shock, twenty2shock);
xlabel('Time');
ylabel('Amplitude');
title('.22 Revolver')
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

%% plot fingerprint t domain

figure
plot(index_t, x_39, 'b');
title('Walter PPQ 9mm Pistol')
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

figure
plot(index_t, q_30, 'r');
title('Ruger 10/22 .22 Rifle')
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

figure
plot(index_t, f_46, 'k');
title('Bersa .380 Pistol');
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)


figure
plot(index_t, x_39, 'b');
hold on
plot(index_t, q_30, 'r');
hold on
plot(index_t, f_46, 'k');
xlim([0 2048])
legend('Walter PPQ 9mm Pistol', 'Ruger 10/22 .22 Rifle', 'Bersa .380 Pistol');
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

figure

subplot(3,1,1)
plot(index_t, x_39, 'b');
legend('Walter PPQ 9mm Pistol');
xlim([0 2048])
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

subplot(3,1,2)
plot(index_t, q_30, 'r');
legend('Ruger 10/22 .22 Rifle');
xlim([0 2048])
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

subplot(3,1,3)
plot(index_t, f_46, 'k');
legend('Bersa .380 Pistol');
xlim([0 2048])
set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)

%% plot fingerprint f domain