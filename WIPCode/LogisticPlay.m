clear
load remedial_reading_data
X = [remedial_reading_data(:,1:2)];
Y = remedial_reading_data(:,3);

[S] = subfnLogisticRegressStats(Y,X);