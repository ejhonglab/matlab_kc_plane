function [S] = FnLoadSyncEpisode(pathandfilename, varargin)
%% Select one h5 file:

[pathname,name,ext] = fileparts(pathandfilename);
filename = strcat(name, ext);

%% Load params from XML:
clockRate = 20000000;
sampleRate = LoadSyncXML(pathname);

%% Start loading HDF5:
%pathandfilename = strcat(pathname,filename);
info = h5info(pathandfilename);

%% Parse input:
props = {'start','length','interval'};
data = {[1,1],[1 Inf],[1 1]};

if(~isempty(varargin))
    assert(rem(length(varargin),2)==0 && iscellstr(varargin(1:2:end)), 'Inputs failed to conform to expected string-value pair format, eg> ''start'',1');
    %foundProps = intersect(varargin(1:2:end), props); 
    IdxCell = cellfun(@(x) strcmpi(x,props),varargin(1:2:end),'UniformOutput',false);
    val = double(cell2mat(varargin(2:2:end)))*sampleRate;
    for i=1:length(val)
        data{cell2mat(IdxCell(i))>0} = [1 val(i)];
    end
end

%% Read HDF5:
S = struct;
for j=1:length(info.Groups)
    for k = 1:length(info.Groups(j).Datasets)
        datasetPath = strcat(info.Groups(j).Name,'/',info.Groups(j).Datasets(k).Name);
        datasetName = info.Groups(j).Datasets(k).Name;
        datasetName(isspace(datasetName))='_';   
        datasetValue = h5read(pathandfilename,datasetPath,data{1},data{2},data{3})';
        % load digital line in binary:
        if(strcmp(info.Groups(j).Name,'/DI'))
            datasetValue(datasetValue>0) = 1;
        end
        % create time variable out of gCtr, 
        % account for 20MHz sample rate:
        if(strcmp(info.Groups(j).Name,'/Global'))
            datasetValue = double(datasetValue)./clockRate;
            datasetName = 'time';
        end
        %assignStr = UniqueName(datasetName);
        %assignin('base',assignStr,datasetValue);
        S.(datasetName) = datasetValue;
    end
end


end




function outStr = UniqueName(str)
%% Generate unique name for variable to be exported.

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
vars = evalin('base','who');
index = 1;
unique = false;
cmpStr = str;

while (~unique)
    ret = cellfun(cellfind(cmpStr),vars);
    if(~any(ret))
        outStr = cmpStr;
        unique = true;
    else
        cmpStr = strcat(str,num2str(index,'%03d'));
        index=index+1;
    end
end

end
