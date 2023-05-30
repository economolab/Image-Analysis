
function stack_array = loadTifStack(parent, fn)

Nfiles = numel(fn);
stack = cell(Nfiles, 1);

% progress bar
k = 1;
files_frac = 0;
update_str = strcat('Loading z-stacks (', num2str(k), '/', num2str(Nfiles), ')');
f = waitbar(files_frac,update_str);
figure(f)

for k = 1:Nfiles
    
    % progress bar
    files_frac = k/Nfiles;
    update_str = strcat('Loading z-stacks (', num2str(k), '/', num2str(Nfiles), ')');
    waitbar(files_frac,f,update_str);

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

stack_array = zeros(Nfiles, size(stack{1}, 2), size(stack{1}, 3), size(stack{1}, 4));
for i = 1:Nfiles

    % progress bar
    files_frac = i/Nfiles;
    update_str = strcat('Squishing z-stacks together (', num2str(i), '/', num2str(Nfiles), ')');
    waitbar(files_frac,f,update_str);

    stack_array(i,:,:,:) = stack{i};
end



close(f)



