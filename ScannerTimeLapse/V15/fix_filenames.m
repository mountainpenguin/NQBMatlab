function fix_filenames(DirName)
%% fix_filenames(DirName)
% -----------------------------------------------------------------------
% Purpose: Fix default filenames from scanner on CMBI2
%
% Description: Renames files from the form 'Name_0_YYYYMMDD_HHMM.tif' to 
%           'P1_00000.tif' where 00000 is the number of minutes since the
%           first image was captured. Also removes the alpha channel from 
%           all images
%
% Arguments: DirName (optional) - Name of directory, defaults to the 
%                                 current working directory.
%
% Input files: 'Name_0_YYYYMMDD_HHMM.tif' - the pictures
% Output files: 'P1_00000.tif' - fixed pictures
% -----------------------------------------------------------------------

if nargin == 0
    %% if no DirName given, use current working dir
    DirName = pwd;
end


disp('Creating backup of /Pictures --> /Pictures.backup');
%% create backup directory for images
indir = fullfile(DirName, 'Pictures');
outdir = fullfile(DirName, 'Pictures.backup');
[status_cp, msg_cp] = copyfile(indir, outdir);
if ~status_cp
    error(msg_cp);
end

disp('Renaming files...');

%% list files
files = dir(fullfile(DirName, 'Pictures', '*.tif'));
t0 = get_time(files(1).name);
for file = files'
    t1 = get_time(file.name);
    tdiff = etime(t1, t0) / 60;
    tdiff = round(tdiff);
    outname = sprintf('P1_%05d.tif', tdiff);
    disp([' ', file.name, ' --> ', outname]);

    im = imread(fullfile(DirName, 'Pictures', file.name));
    im = im(:, :, 1:3);
    inpath = fullfile(DirName, 'Pictures', file.name);
    outpath = fullfile(DirName, 'Pictures', outname);
    %movefile(inpath, outpath);
    delete(inpath)
    imwrite(im, outpath)
end
disp('Done');
end

function t = get_time(name)
f0 = strsplit(name, '_');
f0 = strjoin(f0(3:4), '');
f0 = strsplit(f0, '.');
f0 = char(f0(1));
dstr = sprintf('%s-%s-%s %s:%s', f0(1:4), f0(5:6), f0(7:8), f0(9:10), f0(11:12));
t = datevec(dstr, 'yyyy-mm-dd HH:MM');
end
