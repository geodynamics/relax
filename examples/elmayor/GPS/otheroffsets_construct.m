flocations=fopen('../gps/gps_within200km.dat');
location_data=textscan(flocations,'%s %f %f %f','commentstyle','#');
fclose(flocations);
stations = struct('staNam',[],'east',[],'north',[]);
stations.names = location_data{1};

elmayor_time = 2010.2589;
otheroffsets = cell(length(stations.names),1);

otheroffsets{strcmp(stations.names,'p497')} = [elmayor_time + 0.1945;];
otheroffsets{strcmp(stations.names,'p498')} = [elmayor_time + 0.1945; 2012.16; 2012.27; 2012.55; 2012.75;];
otheroffsets{strcmp(stations.names,'p503')} = [elmayor_time + 0.1945;];
otheroffsets{strcmp(stations.names,'p744')} = [elmayor_time + 0.1945;];
otheroffsets{strcmp(stations.names,'usgc')} = [elmayor_time + 0.1945;];

otheroffsets{strcmp(stations.names,'azry')} = [elmayor_time + 0.2574; 2014.0];
otheroffsets{strcmp(stations.names,'mvfd')} = [elmayor_time + 0.2574;];
otheroffsets{strcmp(stations.names,'pin1')} = [elmayor_time + 0.2574;];
otheroffsets{strcmp(stations.names,'pin2')} = [elmayor_time + 0.2574;];
otheroffsets{strcmp(stations.names,'p479')} = [elmayor_time + 0.2574;];
otheroffsets{strcmp(stations.names,'p486')} = [elmayor_time + 0.2574;];
otheroffsets{strcmp(stations.names,'p505')} = [elmayor_time + 0.2574;];
otheroffsets{strcmp(stations.names,'p741')} = [elmayor_time + 0.2574;];

otheroffsets{strcmp(stations.names,'p500')} = [elmayor_time - 0.2986];
otheroffsets{strcmp(stations.names,'iid2')} = [elmayor_time - 0.2986];

otheroffsets{strcmp(stations.names,'bill')} = [elmayor_time + 2.395; elmayor_time + 2.5;];
otheroffsets{strcmp(stations.names,'bomg')} = [elmayor_time - 0.58; elmayor_time - 0.2; elmayor_time + 0.72;];
%otheroffsets{strcmp(stations.names,'dvlw')} = [elmayor_time + 0.68;];
otheroffsets{strcmp(stations.names,'cact')} = [2014.0; 2014.87];
otheroffsets{strcmp(stations.names,'cotd')} = [2014.0;];
otheroffsets{strcmp(stations.names,'monp')} = [elmayor_time + 2.47;];
otheroffsets{strcmp(stations.names,'nsss')} = [elmayor_time - 6.6;];
otheroffsets{strcmp(stations.names,'p473')} = [elmayor_time - 6.6;];
otheroffsets{strcmp(stations.names,'p491')} = [elmayor_time + 1.89;];
otheroffsets{strcmp(stations.names,'p492')} = [elmayor_time - 0.33;];
otheroffsets{strcmp(stations.names,'p493')} = [elmayor_time + 2.395;];
otheroffsets{strcmp(stations.names,'p494')} = [elmayor_time + 0.4;];
otheroffsets{strcmp(stations.names,'p495')} = [elmayor_time + 0.7;];
otheroffsets{strcmp(stations.names,'p500')} = [elmayor_time + 2.395;];
otheroffsets{strcmp(stations.names,'p507')} = [elmayor_time - 2.0; elmayor_time - 2.25; elmayor_time + 1.9; elmayor_time + 2.395; elmayor_time + 2.55;];
otheroffsets{strcmp(stations.names,'p508')} = [elmayor_time + 2.395];
otheroffsets{strcmp(stations.names,'p741')} = [elmayor_time - 2.45; elmayor_time - 2.1];
otheroffsets{strcmp(stations.names,'p742')} = [elmayor_time + 0.32];
otheroffsets{strcmp(stations.names,'p744')} = [elmayor_time - 1.95; elmayor_time + 2.55;];
otheroffsets{strcmp(stations.names,'pmob')} = [elmayor_time - 5.25];
otheroffsets{strcmp(stations.names,'psap')} = [2014.0;];
otheroffsets{strcmp(stations.names,'sio3')} = [elmayor_time - 5.2];
otheroffsets{strcmp(stations.names,'sio5')} = [elmayor_time + 0.72];
otheroffsets{strcmp(stations.names,'tmap')} = [elmayor_time - 3.6; elmayor_time - 1.9; elmayor_time + 2.3];
otheroffsets{strcmp(stations.names,'widc')} = [2014.0];
otheroffsets{strcmp(stations.names,'wwmt')} = [2014.0; elmayor_time - 4.9];

save('otheroffsets','otheroffsets')