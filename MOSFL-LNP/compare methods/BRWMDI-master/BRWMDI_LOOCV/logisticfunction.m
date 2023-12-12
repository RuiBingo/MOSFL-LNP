function [finds] = logisticfunction(s,c,d)
[row,column]=size(s);
finds=zeros(row,column);
for i=1:row
    for j=1:column
        finds(i,j)=1/(1+exp(c*s(i,j)+d));
    end
end

end
