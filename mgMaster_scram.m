clear all
doSave = 1;

rdmPath = '/data/hammer/space2/mvpaGoals/data/golsaRSA/';
rdmFile = '11-Oct-2017.mat';
load([rdmPath rdmFile]);

clear rdm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Representation Distance Matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    
    %%% Setup %%%
    
    %Possible Regressors: 1=Start, 2=Prompt, 3=fb
    USED_REGS = [1,2,3]; 
    
   
    rdmInfo.trialList = trialList;
    rdmInfo.trialTypes = trialTypes;
    rdmInfo.allTimes = patternTimes;
    rdmInfo.usedRegs = USED_REGS;
    rdmInfo.date = date;
    
    nComponents = size(patternData, 2);
    nTrials = size(trialList, 1);
    nTimes = size(patternTimes, 1);
    nPatterns = nTrials * nTimes;

    
    mgStates = {'cows', 'house', 'scrow', 'stump'};
    mgStates2 = {'Cows', 'House', 'Scrow', 'Stump'};
    mgRegs = {'Start', 'Prompt', 'fb', 'Start', 'Prompt', 'fb'}; 
    
    %get vector telling us which time points to ignore because they are the
    %second-step timepoints in single-step trials
    timeVec = [0, 0, 0, 1 ,1 ,1]; %vector of length nTimes indicating whether each time refers to a second-trial time point
    if length(timeVec) ~= nTimes
        error('timeVec wrong length')
    end
    colTimeTypes = repmat(timeVec, 1, nTrials); 
    trialTypeTemp = repmat(trialTypes, 1, nTimes);
    colTrialTypes = reshape(trialTypeTemp', 1, nTimes * nTrials);

    colRemove = (colTimeTypes==1) & (colTrialTypes==1);
    
    %Add columns from unused regressors to be removed
    startCols = repmat([1, 0, 0], 1, nPatterns/3);
    promptCols = repmat([0, 1, 0], 1, nPatterns/3);
    fbCols = repmat([0, 0, 1], 1, nPatterns/3);
    
    usedRegCols = startCols * any(USED_REGS==1) + promptCols * any(USED_REGS==2) + fbCols * any(USED_REGS==3);
    
    colRemove(~usedRegCols) = 1;
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% create ordered list of beta names corresponding to patterns %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numBetas =  nTrials * nTimes; 
    betaNames = {}; %list of filenames
    betaTrials = {};
    sameBetas=[];
    for beta_i = 1:numBetas


        trialPath = trialList(floor((beta_i - .5)/nTimes)+1,:); %e.g. [2,4,3]
        timePoint = mod(beta_i, nTimes);
        
        if timePoint == 0
            timePoint = nTimes;
        end
        
        if trialPath(2)==trialPath(3)
            finalGoal = 'None';
        else
            finalGoal = mgStates2{trialPath(3)};
        end
        
        if timePoint <= 3
            transition = [mgStates{trialPath(1)} mgStates2{trialPath(2)} finalGoal];
            
        else
            transition = [mgStates{trialPath(2)} mgStates2{trialPath(3)} finalGoal];
        end
        
        if ~any(timePoint == [3, 6]) %if not feedback
            betaName = [transition mgRegs{timePoint}]; %eg cowsHouseStumpStart
        else
            betaName = [mgRegs{timePoint} transition]; %eg fbcowsHouseStump
        end
        
        %if current beta name is already in betaNames and won't be removed
        %at the end due to colRemove, it's a duplicate so store it in list
        %of 'sameBetas'. The patterns corresponding to the same betas will
        %be averaged
        if any(cellfun('isempty',strfind(betaNames,betaName))==0) && ~any(beta_i == find(colRemove)) %&& timePoint==4
%             disp('==============================')
%             beta_i
%             betaName
%             trialPath
%             
%             oldIdx = find(cellfun('isempty',strfind(betaNames,betaName))==0)
%             disp('Old One')
%             betaTrials{oldIdx}
%             betaNames{oldIdx}
%             
%             disp('==============================')

            oldIdx = find(cellfun('isempty',strfind(betaNames,betaName))==0);
            sameBetas = vertcat(sameBetas, [oldIdx, beta_i]);
        end
        
        betaNames{beta_i} = betaName;
        betaTrials{beta_i} = trialPath;            
    end
    
    
    %remove columns corresponding to duplicate betas
    % e.g. the second half of [1,2,4] and [4,2,4]
    if ~isempty(sameBetas)
        colRemove(sameBetas(:,2)') = 1;
    end
    
    betaNames(colRemove) = [];
    rdmInfo.betaNames = betaNames;

    
    
    %for debugging; make sure all betaNames are unique
    [U, ic, iu] = unique(betaNames);
    count = histc(iu, 1:numel(ic));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% reshape pattern data for RDM %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for comp_i = 1:nComponents
        
        %average across repeats
        compDataMean = mean(patternData(comp_i).data, 4);  %patternData = [timepoints, units, trials, reps]
        
        %reshape to have one column for every timepoint, trialtype pair
        %columns representing timepoints from the same trial are next to
        %each other
        
        compData = nan([size(compDataMean,2), size(compDataMean,1) * size(compDataMean,3)]); %[units, timepoint x trial]
        for trial_i=1:size(compDataMean,3)
            colStart = (trial_i - 1) * nTimes + 1;
            colEnd = trial_i * nTimes;
            compData(:, colStart:colEnd) = compDataMean(:,:,trial_i)';
        end
        
        %%%%%%%%%%%%%%%%%% 
        %%% Clean Data %%%
        %%%%%%%%%%%%%%%%%%
        
        %Average data across trials with the same betas (see analysis
        %notes 10/10/17)
        for sameBeta_i = 1:size(sameBetas, 1)
            betaPair = sameBetas(sameBeta_i, :);
            
            compData(:, betaPair(1)) = (compData(:, betaPair(1)) + compData(:, betaPair(1)))/2;
            
        end
              
        
        %Remove second-trial times from one-step trials
        compDataClean=compData;
        compDataClean(:,colRemove)=[];
      

        %%%%%%%%%%%%%%%%%%%%%%%%% 
        %%% Compute Final RDM %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        compRDM = 1 - corr(compDataClean);
        rdm.(layers{comp_i}.name) = compRDM; %corrected to be 1-corr, 10/9/17
%         
    end
    
    
    %save component similarity matrix
    
    
    
    outputRoot = '/data/hammer/space2/mvpaGoals/data/golsaRSA/';
    outputFile = [outputRoot date '.mat'];
    
    
            
    if doSave
        if exist(outputFile, 'file')
            resp = input(['File ' outputFile ' already exists. Overwrite? (Y/N)'], 's');
            if any(resp ==  ['y', 'Y'])
                save(outputFile, 'rdmInfo', 'rdm', 'patternData') 
            else
                disp('File not written.')
            end
        end     
    end
    
end

%double check patterns and rdms
if 0
    
    
    
    
 
    
    %use inspectPatterns.m to look at pattern data
%     allLayers = fieldnames(rdm);
%     checkLayer='state';
%     layerNum = find(cellfun('isempty',strfind(allLayers,checkLayer))==0);
%     nTimes = size(patternTimes, 1);
%     
%     checkVals = corr(checkData);
%     timePoint = 1;   
%     checkVals = checkVals(timePoint:nTimes:end, timePoint:nTimes:end);
%     figure(88)
%     imagesc(checkVals)
end
    

