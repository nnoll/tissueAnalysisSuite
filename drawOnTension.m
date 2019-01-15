t = 20;

T = PN{1}.returnTension();
Tmax = prctile(T,95);
Tmin = prctile(T,5);

T = (T-Tmin)/(Tmax-Tmin);
T(T>1) = 1;
T(T<0) = 0;

% Convert tension to color.
cmap = hot(256);
x = linspace(0,1,256);
clear Tcolor
Tcolor(:,1) = interp1(x,cmap(:,1),T);
Tcolor(:,2) = interp1(x,cmap(:,2),T);
Tcolor(:,3) = interp1(x,cmap(:,3),T);

r = PN{1}.computePrimalVerts();

hold on
r1 = zeros(size(PN{1}.d1,2),2);
r2 = zeros(size(PN{1}.d1,2),2);
badEdges = [];
for e = 1:size(PN{1}.d1,2)
    verts = find(PN{1}.d1(:,e)~=0);
    if (length(verts) == 2)
        r1(e,:) = r(verts(1),:);
        r2(e,:) = r(verts(2),:);
    else
        badEdges = [badEdges,e];
    end
end

Tcolor(badEdges,:) = [];
r1(badEdges,:) = [];
r2(badEdges,:) = [];

% img = mem(:,:,t);
rgb = cat(3,mat2gray(mem(:,:,t)),mat2gray(mem(:,:,t)),mat2gray(mem(:,:,t)));
img = zeros(size(mem,1),size(mem,2),3);
img = insertShape(double(img),'Line',[r1,r2],'Color',Tcolor,'Opacity',1);
img = imdilate(img,strel('disk',1));
% rgb(img>0) = 0;
rgb = .5*rgb + img;

close all
imshow(rgb)

