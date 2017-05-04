function [X,Y,max_distance,sequence]=fun_find_max_distance(pressure,hole,speed,distance)
% This function only works for data type 3 (3 rows of order variable)
X.origin=hole;
Y.origin=speed;
P=length(pressure);
m=length(hole);
n=length(speed);
k=1;
[rows,cols]=size(distance);
if (rows >=3)
    error('Wrong Order type of Input')
    break
end
sequence(1:rows,:)=distance;
for p=1:P
    for i=1:m
        for j=1:n
            sequence(5,k)=k;
            max_distance.origin(i,j,p)=distance(4,k);
            k=k+1;
        end
    end
end

pressure2=sort(pressure);
X.seq=sort(hole);
Y.seq=sort(speed);
k=1;
index=sequence(rows+1,:);

for p=1:P
    fp= find(sequence(1,index)==pressure2(p));
    for i=1:m
        fh= find(sequence(2,index(fp))==X.seq(i));
        for j=1:n
            fs= find(sequence(3,index(fp(fh)))==Y.seq(j));
            kk=index(fp(fh(fs)));
            max_distance.seq(i,j,p)=sequence(4,kk);
            sequence(rows+2,k)=kk;
            k=k+1;
            clear fs
        end
        clear fh
    end
    clear fp
end

