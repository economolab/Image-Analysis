% get file names for all files of a certain type in a directory
% returns a cell array of file names
% takes parent (string): directory of interest, and filetype (string): file
% types to be extracted, e.g. tif, jpg, etc. 

function fns = getFileNames(parent, filetype)

fns = {};

all_files = extractfield(dir(parent), 'name');
for i = 1:length(all_files)

    if contains(all_files{i}, strcat('.', filetype))
        fns{end+1} = all_files{i}; %#ok
    end

end

