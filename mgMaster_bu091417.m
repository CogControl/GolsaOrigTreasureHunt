ssstrialList=vertcat([1,2,4], [1,2,2], [1,3,4], [1,3,3],...
    [2,1,3], [2,1,1], [2,4,3], [2,4,4],...
    [3,1,2], [3,1,1], [3,4,2], [3,4,4],...
    [4,3,1], [4,3,3], [4,2,1],[4,2,2]);
trialList = trialList(1:4,:);
numReps = 10; %number of times to repeat each trial type

results = repmat(horzcat(trialList, zeros(size(trialList,1))), [1,1,numReps]); %ensure correct movement (for db)


patternTimes = [
    5,7; %loading part 1
    17,19; %loading part 2
    20, 22; %maintaining both
    24, 26; %acting
    26, 28]; %outcome

%initialize network
%clear n
%clear e
mgNet;
% mgRunTrial(n,e,4,2,4,10,10,0);
%mgRunTrial(n,e,1,2,4)

%initialize results struct
layers = n.layerList;

for n_i = 1:length(layers)
    patternData(n_i).name = layers{n_i}.name;
    patternData(n_i).data = nan(size(patternTimes,1), layers{n_i}.numUnits,...
    size(trialList,1), numReps);
end
    

for trial_i = 1:size(trialList, 1)
    trial = trialList(trial_i, :);
    trial_i
    
    for rep_i = 1:numReps
        endLoc = mgRunTrial(n, e, trial(1), trial(2), trial(3));
        results(trial_i, 4, rep_i) = endLoc;

        if endLoc ~= trial(2)
            trial_i
            trial
            error('wrong state reached')
        end

        for layer_i = 1:numel(layers)
            layer = layers{layer_i}; 
            newData = n.getPatterns(layer.name, patternTimes);
            patternData(layer_i).data(:, :, trial_i, rep_i) = newData;
        end
    end
    
end

%add results column - 1 if state transition correct, else 0
results(:,5,:) = results(:,2,:) == results(:,4,:);


%calculate variance across trials

variances = zeros(size(layers, 2), size(trialList, 1));

for layer_i = 1:size(layers,2)
    for trial_i = 1:size(trialList, 1)
        curData = patternData(layer_i).data(:,:, trial_i,:);
        
        variances(layer_i,trial_i) = mean(mean(var(curData,0,4),1),2);
    end
end



if 0
    %make  similarity matrices
    nPatterns = size(trialList, 1) * size(patternTimes, 1);
    nComponents = size(patternData, 2);
    simMatrices = zeros(nPatterns, nPatterns, nComponents);
    for comp_i = 1:size(patternData,2)
        compData = mean(patternData(comp_i).data, 4);
        compData2 = reshape(compData, [size(compData,1), size(compData,2) * size(compData,3)]);
        corrData = corr(compData2);
        figure(44)
        hold on
        
        imagesc(corrData)
        title(patternData(comp_i).name)
        hold off
        
        pause()
        close(44)
        %corr
    end
end

