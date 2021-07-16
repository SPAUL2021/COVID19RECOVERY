
nd = 14;

for i = 1:nd
    tau(i) = i; % (1/3)+(i-1);
end   

ns = 29;

for i = 1:ns
    sigma(i) = 14+(i-1); % + (1/3);
end  

data = load('dts.txt');

for i = 1:nd
    for j = 1:ns
        rptm(j,i) = data(ns*(i-1)+j);
    end
end    


   init = 0;
for i = 1:nd
    for j = 1:ns
        test = rptm(j,i);
        if (test >= 10)
            n = round(test/10);
            for m = 1:n
                x(init+m) = i;
                y(init+m) = j+13;
            end    
                init = init+n;
        end   
    end
end

s1 = size(x)
s2 = size(y)
init9 = init

inc = zeros(init,1);
rec = zeros(init,1);

for i = 1:init
    inc(i) = x(i);
    rec(i) = y(i);
end

 A = [inc; rec];

 fileID = fopen('bivariD.txt','w')
 fprintf(fileID,'%12s\n','size');
 fprintf(fileID,'%d\n',init);
 fprintf(fileID,'%32s\n','incubation(Column_1)');
 fprintf(fileID,'%d\n',inc);
 fprintf(fileID,'%32s\n','recovery(Column_2)');
 fprintf(fileID,'%d\n',rec);


%------------------------------------------------------------------