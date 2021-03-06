function [tot_mhour, tot_probe, tot_temperature, tot_p_type] = ReadThermometer(FullFileName)
% [tot_mhour, tot_probe, tot_temperature, tot_p_type] =
%       ReadThermometer(FullFileName)
% -------------------------------------------------------------------------
% Purpose : reading a file generated by the thermometer
% Description : The file looks like this:
%       HEADERS
%       ...
%       DATE: 01.01
%       00:00:00 probe temp type
%                probe temp type
%                   ...
%       00:00:30 probe temp type
%                   ....
%
% Arguments : FullFileName - The name of the thermometer file
% Returns : All parameters are cell arrays the size of different 
%       experiments made.
%       tot_mhour - each cell contains a vector of time of measurment
%       tot_probe - each cell contains the list of probes
%       tot_temperature - each cell contains a matrix of temperature per
%       probe
%       tot_p_type - each cell contains the list of types of probes
% -------------------------------------------------------------------------
% Irit Levin. 10.01.08

% showing the "open file" ui, if it wasn't specified
if nargin == 0
    [FileName, FilePath] = uigetfile('*.txt');
    if isequal(FileName,0)
        return;
    end
    FullFileName = [FilePath, FileName];
end

% reading the data
% -----------------

% after the first time 'DATA:' appears, the data begins
fid = fopen(FullFileName);
str = fgetl(fid);
k   = [];
while (~size(k,1) && ~feof(fid))
    k = strfind(str, 'DATE');
    str = fgetl(fid);
end

hr          = '00:00:00';
mhour       = [];
probe       = [];
temperature = [];
p_type      = [];
NDataSets   = 1;
NProbe      = 1;
NLine       = 0;
% after the headers - reading the data
while ~feof(fid)
    k = strfind(str, 'DATE');
    if size(k,1)
        tot_mhour{NDataSets} = mhour;
        tot_probe{NDataSets} = probe;
        tot_temperature{NDataSets} = temperature;
        tot_p_type{NDataSets} = p_type;
        NDataSets = NDataSets+1;
        mhour       = [];
        probe       = [];
        temperature = [];
        p_type      = [];
        NLine       = 0;
    else
        if isspace(str(1))
            % format is: ____ probe temperature type
            s  = textscan(str, '%2s: %f �C %s');
            %mhour{NLine, NProbe}    = hr;
            probe{NProbe}          = char(s{1});
            temperature(NLine, NProbe) = s{2};
            p_type{NProbe}         = s{3};
            NProbe = NProbe + 1;
        else
            % format is: time probe temperature type
            NLine = NLine+1;
            s = textscan(str, '%8s %2s: %f �C %s');
            NProbe  = 1;
            hr      = char(s{1});
            mhour{NLine}           = hr;
            probe{NProbe}          = char(s{2});
            temperature(NLine, NProbe) = s{3};
            p_type{NProbe}         = s{4};
            NProbe = NProbe + 1;
        end
    end
    str = fgetl(fid);
end
tot_mhour{NDataSets} = mhour;
tot_probe{NDataSets} = probe;
tot_temperature{NDataSets} = temperature;
tot_p_type{NDataSets} = p_type;

fclose(fid);