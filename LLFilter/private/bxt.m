function b = bxt(x,g,Jg)
m=size(Jg,2);
for i=1:m
   b{i}=g{i}-Jg{i}*x;
end