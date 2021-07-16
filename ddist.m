


data = load('dsig.txt');

n1 = size(data);
n = n1(1)

mdata = min(data);

rdata = (1/mdata)*data;

N = int16(rdata);


fileID = fopen('ddist.txt','w')

sum = 0;

for j = 1:n
    sum = sum + N(j);
end

sum1 = sum

d = zeros(sum,1);
    ns = 0;
for j = 1:29
    for i = 1:N(j)
        d(ns+i) = 13 + j;
    end
    ns = ns+N(j);
end

d1 = size(d)

fprintf(fileID,'%d\n',d);
fclose(fileID);
