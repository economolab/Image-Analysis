function [clips, maxvals, ix] = getClips(im, c, clipsz)
ypix = size(im, 1);
xpix = size(im, 2);
Nchan = size(im, 3);
Nplanes = size(im, 4);

if mod(clipsz,2)==0
    clipsz = clipsz+1;
end
cliphf = floor(clipsz/2);

clips = zeros(clipsz, clipsz, c.N, Nchan, Nplanes);

ix.minx = zeros(c.N, 1);
ix.maxx = zeros(c.N, 1);
ix.miny = zeros(c.N, 1);
ix.maxy = zeros(c.N, 1);

for i = 1:c.N
    x1 = c.x(i)-cliphf;
    x2 = c.x(i)+cliphf;
    y1 = c.y(i)-cliphf;
    y2 = c.y(i)+cliphf;

%     z = c.z(i);

    cx1 = 1;
    cy1 = 1;
    cx2 = clipsz;
    cy2 = clipsz;

    if x1<1
        x1 = 1;
        cx1 = abs(c.x(i)-cliphf)+2;
    end
    if y1<1 
        y1 = 1;
        cy1 = abs(c.y(i)-cliphf)+2;
    end   
    if x2>xpix
        x2 = xpix;
        cx2 = 1+x2-x1;
    end
    if y2>ypix
        y2 = ypix;
        cy2 = 1+y2-y1;
    end

    for j = 1:Nchan
        
        ix.minx(i) = x1;
        ix.maxx(i) = x2;
        ix.miny(i) = y1;
        ix.maxy(i) = y2;

        clips(cy1:cy2, cx1:cx2, i, j, :) = im(y1:y2, x1:x2, j, :);
        for k = 1:Nplanes
            clips(:,:,i,j, k) = imgaussfilt(clips(:,:,i,j,k),1);

        end

    end
end

maxvals = zeros(Nchan, 2);
for i = 1:Nchan
    tmp = clips(:,:,:,i, :);
    maxvals(i, 1) = prctile(tmp(:), 25);
    maxvals(i, 2) = prctile(tmp(:), 99.9);
end


