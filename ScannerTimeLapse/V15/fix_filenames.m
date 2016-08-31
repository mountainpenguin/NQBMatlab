function fix_filenames(DirName, timestep)
%% fix_filenames(DirName, timestep)
% -----------------------------------------------------------------------
% Purpose: Fix default filenames from scanner on CMBI2
%
% Description: Renames files from the form 'Name_0_YYYYMMDD_HHMM.tif' to 
%           'P1_00000.tif' where 00000 is the number of minutes since the
%           first image was captured.
%
% Arguments: DirName (optional) - Name of directory, defaults to the 
%                                 current working directory.
%            timestep (optional) - Time between images, defaults to 15 
%                                  minutes.
%
% Input files: 'Name_0_YYYYMMDD_HHMM.tif' - the pictures
% Output files: 'P1_00000.tif' - fixed pictures
% -----------------------------------------------------------------------

if nargin == 0
    %% if no DirName given, use current working dir
    DirName = pwd;
end

if nargin < 2
    %% if no timestep given, use 15
    timestep = 15;
end

%% create backup directory for images
[success_mk, msg_mk, msgid_mk] = mkdir(DirName, 'Pictures.backup');
if ~success_mk
    error(msgid_mk, msg_mk);
end

disp('-----------------------------------------------------------------');
disp([datestr(now), ' ', DirName]);
disp('-----------------------------------------------------------------');
disp('Renaming files...');

%% list files
files = dir(fullfile(DirName, 'Pictures', '*.tif'));
t0 = datevec(files(1).date, 'dd-mmm-yyyy HH:MM:SS');
for file = files'
    t1 = datevec(file.date, 'dd-mmm-yyyy HH:MM:SS');
    tdiff = etime(t1, t0) / 60;
    tdiff = round(tdiff);
    outname = sprintf('P1_%05d.tif', tdiff);
    %disp([' ', file.name, ' --> ', outname]);
    inpath = fullfile(DirName, 'Pictures', file.name);
    outpath = fullfile(DirName, 'Pictures', outname);
    movefile(inpath, outpath);
end
disp('Done');
end
