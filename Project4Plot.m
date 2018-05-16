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
Neurotoxin = logical(GrainDeerData(:,7));

neuroMonth = month(Neurotoxin);
neuroDeer = GrainDeer(Neurotoxin);

figure
plot(month,GrainHeight,'b-')
hold on
plot(month,GrainDeer,'r-')
plot(neuroMonth,neuroDeer,'k.')
xlabel('Month')
%legend('Temperature, Celsius', 'Precipitation, cm', 'Grain Height, cm', 'Grain Deer Population', 'Neurotoxin')
legend('Grain Height, in', 'Grain Deer Population', 'Neurotoxin Event')
legend('location', 'eastoutside')
axis([-inf,inf,0,inf])
hold off

figure
plot(month(1:60),Temperature(1:60))
hold on
plot(month(1:60),Precipitation(1:60))
xlabel('Month')
legend('Temperature, Celsius', 'Precipitation, cm')
legend('location', 'eastoutside')
hold off