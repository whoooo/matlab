index = zeros(10,3);
index(:,1) = linspace(1,10,10);
index(:,2) = linspace(20,30,10);
index(:,3) = linspace(50,60,10);

data = zeros(10,3);
data(:,1) = 1;
data(:,2) = 2;
data(:,3) = 1;

plot(index,data);