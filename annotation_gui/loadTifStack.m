
function stack = loadTifStack(parent, fn)

Nfiles = numel(fn);
stack = cell(Nfiles, 1);


for k = 1:Nfiles
    fullfn = fullfile(parent, fn{k});
    info = imfinfo(fullfn);

    Nrow = info(1).Height;
    Ncol = info(1).Width;

%     [Nchans, Nplanes] = parseInfo(info(1));
    Nchans = numel(info(1).MaxSampleValue);
    Nplanes = numel(info);

    stack{k} = zeros(Nchans, Nrow, Ncol, Nplanes);

    cnt = 0;
    for i = 1:Nplanes
        for j = 1:Nchans
            cnt = cnt+1;
            stack{k}(j,:,:,i) = imread(fullfn, 'tif', 'Info', info, 'Index', cnt);
%             stack{k}(:,:,:,i) = permute(imread(fullfn, 'tif', 'Info', info, 'Index', cnt), [3, 1, 2]);

        end
    end

end
stack = cell2mat(stack);
stack = permute(stack, [2 3 1 4]); 
