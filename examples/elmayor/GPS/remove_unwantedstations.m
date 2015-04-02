fstations=fopen('gps_locations_km.dat');
stations=textscan(fstations,'%s %f %f %f','commentstyle','#');
fclose(fstations);

inside_names = stations{1}(sqrt(stations{2}.^2 + stations{3}.^2)<200);
inside_east = stations{2}(sqrt(stations{2}.^2 + stations{3}.^2)<200);
inside_north = stations{3}(sqrt(stations{2}.^2 + stations{3}.^2)<200);

fexport=fopen('gps_within200km.dat','wt');
for k = 1:length(inside_names)
    fprintf(fexport,'%s %f %f \n',inside_names{k},inside_east(k),inside_north(k));
end
fclose(fexport);

outside_names = stations{1}(sqrt(stations{2}.^2 + stations{3}.^2)>=200);
for k = 1:length(outside_names)
    delete(['./put_SOPAC_timeseries_here/' outside_names{k} 'FilterTrend.neu'])
end