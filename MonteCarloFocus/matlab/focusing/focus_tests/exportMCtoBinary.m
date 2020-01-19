

file='data2fluoro.mc';

fileID=fopen(file,'w');

fwrite(fileID,Nphotons,'uint32');
fwrite(fileID,p,'float');
fwrite(fileID,v,'float');
fwrite(fileID,p0,'float');
fwrite(fileID,w,'float');

fclose(fileID);