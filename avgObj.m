fid=fopen('OBJFUC.DAT','r');   
linn=fgets(fid);
linn=fgets(fid);
linn=fgets(fid);

obj = []
while ~feof(fid)
    linn=fgets(fid);
    temp = sscanf(linn,'%f');
    obj = [obj; temp(2)];
end

fclose(fid);