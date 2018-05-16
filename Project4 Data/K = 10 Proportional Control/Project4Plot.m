load('GrainDeerData.txt')

month = GrainDeerData(:,1);
numYears = GrainDeerData(end,2) - GrainDeerData(1,2);
year = [1:12:12*numYears + 1];
year = repmat(year, 12, 1);
year = reshape(year, [], 1);
month = month + year;
Temperature = GrainDeerData(:,3);
Precipitation = GrainDeerData(:,4);
GrainHeight = GrainDeerData(:,5)./2.54;
GrainDeer = GrainDeerData(:,6);
Neurotoxin = GrainDeerData(:,7);

figure
plot(month,GrainHeight,'b-')
hold on
plot(month,GrainDeer,'r-')
plot(month,Neurotoxin,'k-')
xlabel('Month')
legend('Cake Grain Height, in', 'Frail Grain Deer Population', 'Deadly Neurotoxin Units')
legend('location', 'eastoutside')
axis([-inf,inf,0,142])
hold off

figure
plot(month(1:60),Temperature(1:60))
hold on
plot(month(1:60),Precipitation(1:60))
xlabel('Month')
legend('Temperature, Celsius', 'Precipitation, cm')
legend('location', 'eastoutside')
hold off