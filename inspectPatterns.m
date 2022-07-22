%inspect patterns used in RSA

if 1
    clear all

    %load saved rdm data
    rdmPath = '/data/hammer/space2/mvpaGoals/data/golsaRSA/';
    rdmFile = '19-Oct-2017.mat';
    load([rdmPath rdmFile]);
    layerNames = fieldnames(rdm);

end

%assumes rdm and patternData use same order of field names

layerName = 'state';
trialNum = 24;
timePoint = 4;

layerIdx = find(ismember(layerNames,layerName));

if ~strcmp(patternData(layerIdx).name, layerName)
    error('What is happening???')
end

patterns = patternData(layerIdx).data;
patterns = mean(patterns,4);

pattern = patterns(timePoint, :, trialNum);

plotData = repmat(pattern, 10, 1);

figure()
plot(plotData)
disp(rdmInfo.trialList(trialNum, :))


