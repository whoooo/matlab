x = linspace(1,14,14);
n = 2.^x;

fast = 4.*n.*log2(2.*n) + n;
slow = n.^2;

figure
a = semilogy(n,slow,n,fast);
% title('Computations for Fast vs Normal Cross Correlation')
legend('Normal (Time Domain)','Fast (Frequency Domain)');
xlabel('N Points')
ylabel('Computations required')
grid minor

set(gca,'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18)