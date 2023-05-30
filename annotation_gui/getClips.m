function ix = getClips(im, c, clipsz)

ypix = size(im, 2);
xpix = size(im, 3);
Nchan = size(im, 1);
Nplanes = size(im, 4);

if mod(clipsz,2)==0
    clipsz = clipsz+1;
end
cliphf = floor(clipsz/2);

ix.minx = zeros(c.N, 1);
ix.maxx = zeros(c.N, 1);
ix.miny = zeros(c.N, 1);
ix.maxy = zeros(c.N, 1);

for i = 1:c.N

    x1 = c.x(i)-cliphf;
    x2 = c.x(i)+cliphf;
    y1 = c.y(i)-cliphf;
    y2 = c.y(i)+cliphf;

    if x1<1
        x1 = 1;
    end
    if y1<1 
        y1 = 1;
    end   
    if x2>xpix
        x2 = xpix;
    end
    if y2>ypix
        y2 = ypix;
    end

    for j = 1:Nchan
        
        ix.minx(i) = x1;
        ix.maxx(i) = x2;
        ix.miny(i) = y1;
        ix.maxy(i) = y2;

    end

end



