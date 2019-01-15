imshow(mem)
hold all
plot.curveSkel(cStruct,'r')

[x,y] = ginput;
x = round(x);
y = round(y);
if (~exist('hairCells'))
    hairCells = unique(L(y + size(L,1)*(x-1)));
else
    hairCells = [hairCells;unique(L(y + size(L,1)*(x-1)))];
end