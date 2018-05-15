load('GrainDeerData.txt')

month = GrainDeerData(:,1);
year = [1:12:61];
year = repmat(year, 12, 1);
year = reshape(year, [], 1);
month = month + year;
Temperature = GrainDeerData(:,3);
Precipitation = GrainDeerData(:,4);
GrainHeight = GrainDeerData(:,5);
GrainDeer = GrainDeerData(:,6);

figure
plot(month,Temperature)
hold on
plot(month,Precipitation)
plot(month,GrainHeight)
plot(month,GrainDeer)
xlabel('Month')
legend('Temperature, Celsius', 'Precipitation, cm', 'Grain Height, cm', 'Grain Deer Population')
axis([-inf,inf,0,80])
hold off