classdef network < handle

% Network object containing all layers, projections, nodes, and handling update
%     
% 
% Properties
%    layerList - cell array of layers (objects) of model neurons
%    projectionList - cell array of projections (objects) between layers of model 
%        neurons (each projection is an object)
%    nodeList - cell array of nodes (objects) which have an output but
%        aren't modeled neurally. 
%   displayList - cell array of objects whose contents should be displayed
%         in "real time" as the model simulation advances
%   accessDict - map (dictionary) from object names to the objects
% 
% MethodsweightIM

%      network - class constructor, sets up empty network
%      addLayer, addNode - adds component to network using specified parameters
%         (calls constructor for relevant class and adds component to name to
%          relevant list).

%       so - selects object by name 
%       
%       update - Main network update function. Outer loop cycles through
%             layers, inner loops cycle through projections to that layer,
%             nodes connected to that layer, and inhibitors connected to
%             that layer. Finally, nodes are cyclically updated as well.
%             
%       connect - connects a layer or node to a layer or projection, 
%       creating a new projection
%       
%       plot - plots timecourse of activity for all specified layer objects 
%           over the specified time
          
%     Network methods serve to add network components, connect network
%     components (via a projection) or update the network. Details on
%     arguments for all methods can be found below.


    
   
    properties
        layerList = {}
        projectionList = {}
        nodeList = {}
        displayList = {{},{}}
        accessDict = containers.Map
        plotType = ''
        
    end
    
    methods
        function obj = network()         
        end
 
        function obj = addLayer(obj, name, selfParams) 
           %selfParams = {numUnits, actFunc, timeConst, decayRate, traceType, inhibitors}
           %inhibitors are layers or nodes that strongly inhibit unit activity
           %see layer.m for more detail. 
            try obj.accessDict(name);                
                error('A layer with that name already exists')
            end
            switch numel(selfParams)
                case 4
                    newLayer = layer(name,selfParams{1},selfParams{2},selfParams{3},selfParams{4});
                case 5
                    newLayer = layer(name,selfParams{1},selfParams{2},selfParams{3},selfParams{4},selfParams{5});
                case 6
                    newLayer = layer(name,selfParams{1},selfParams{2},selfParams{3},selfParams{4},selfParams{5},selfParams{6});
            end
            obj.layerList{end+1} = newLayer;
            obj.accessDict(name) = newLayer;
            
        end
        
        %nodes provide non-neural inputs and signals
        function obj = addNode(obj,name,numOutputs,updateFunc,inputs)
            
            try obj.accessDict(name);
                error('A node with that name already exists')
            end
            if nargin>3
                newNode = node(name,numOutputs,updateFunc,inputs);
            else
                newNode = node(name,numOutputs);
            end
            
            obj.nodeList{end+1} = newNode;
            obj.accessDict(name) = newNode;
        end
        
        function obj=addOscNode(obj,name, timeFactor)
            try obj.accessDict(name);
                error('A node with that name already exists')
            end
            newNode = oscNode(name, timeFactor);
            obj.nodeList{end+1} = newNode;
            obj.accessDict(name) = newNode;
        end
            
       
        
        %method to access network components (select object) by name
        %projection and node names assumed to start with p_ and n_,
        %respectively
        function dobj = so(obj,name)
            try
                dobj = obj.accessDict(name);
            catch
                dobj = 0;
                error(['Error: no network component named ' name])
                return
            end   
        end
        
        
   
        function obj = setNode(obj,nodeName,vals)
            obj.so(nodeName).set_vals(vals);
        end
  
        
        function obj = update(obj)
            %Main network update function. Outer loop cycles through
            %layers, inner loops cycle through projections to that layer,
            %nodes connected to that layer, and inhibitors connected to
            %that layer. Finally, nodes are cyclically updated as well.
            
            global t dt
 
            for layer_i = 1:numel(obj.layerList) %cycle through layers
                curLayer = obj.layerList{layer_i};                     
                
                sumExcit = zeros(1,curLayer.numUnits);
                sumInhib = zeros(1,curLayer.numUnits);
                
                projs = curLayer.get_projections();
                
%                 if strcmp(curLayer.name, 'conjOut')
%                     error()
%                 end
                    
                for proj_i = 1:numel(projs) %cycle through projections to calculate input to layer

                    curProj = projs{proj_i};                       
                    numGates = numel(curProj.gates);
                    
                    
                    %Nodes and layer units can act as gates for a particular
                    %projection. Each gate has a condition which it's value
                    %has to meet for the gate to be open (e.g. .5<val<.8).
                    %All gates must be open for signal to flow through the 
                    %projection. Gates open or close based on a certain 
                    %condition. This loop checks to see if any gates are closed.
                    gatesOpen = true;
                    for gate_i = 1:numGates
                        
                        if ischar(curProj.gates{gate_i}{1}) %if the learnGate object is just the name, replace it with the object (backward compatibility issue)
                            curProj.gates{gate_i}{1} = obj.so(curProj.gates{gate_i}{1});        
                        end
                        
                        %curGate = obj.so(curProj.gates{gate_i}{1}); %curGate = {gateName,polarity}
                        curGate =  curProj.gates{gate_i}{1};
                        gateCondition = curProj.gates{gate_i}{2}; %could be singleton -1 or 1, or could be an interval where the gate is open
                        gateVal = curGate.get_vals();
                        gateDt = curGate.get_dt();
                        if length(gateCondition) == 1
                           
                            if curGate.get_vals() == 0
                                gateVal = -1;
                            else
                                gateVal = 1;
                            end
                            
                            if gateVal*gateCondition == -1
                                gatesOpen = false;
                            end
                            
                        else
                            v1 = min(gateCondition(1:2)); v2 = max(gateCondition(1:2));
                            if ~(gateVal >= v1 && gateVal <= v2) %is oscillation in the appropriate phase?
                                gatesOpen = false;
                            elseif length(gateCondition) == 3 %if in appropriate interval, check for dt condition
                                if (gateDt*gateCondition(3))<0 
                                    gatesOpen = false;
                                end
                            end

                        end
                    
                                   
                        
                    end
                    
                    %if any of a projections gates negative (ie the gate is 'up'), no new activity
                    %can pass through
                    if gatesOpen
                        sourceAct = curProj.get_source().get_act();
                        weightedAct = sourceAct*curProj.get_weights();
                    else                 
                        %sourceAct = zeros(1,curLayer.numUnitsgate);
                        sourceAct = curProj.get_source().get_act(); %source act needs to be >0 so that stateSim->adjFuture can learn even when gate off
                        weightedAct = zeros(1,curLayer.numUnits);
                    end
                    
                    
                    %channels only allow input through units corresponding
                    %to the active units of a particular modulator. see
                    %projection for more details
                    if ~isempty(curProj.channels) %channel = {name, threshold}
                        channelVals = curProj.channels{1}.get_act()>curProj.channels{2};
                        weightedAct = weightedAct.*channelVals;
                    end                     
                    
%                     if atTime(2.0) && curProj==obj.so('p_conjDes_qIn')
%                         error('abc')
%                     end
                    
                    sumExcit = sumExcit+max(weightedAct,0);
                    sumInhib = sumInhib+min(weightedAct,0);


                    if ~isempty(curProj.modulators)
                        for mod_i  = 1:numel(curProj.modulators)
                            modulator = curProj.modulators{mod_i};
                            if ischar(modulator)
                                modulator=obj.so(modulator);
                                curProj.modulators{mod_i} = modulator;
                            end
                                
                            curProj.modVals(mod_i,:) = modulator.act;
                        end
                    end

                    %update weights if learning is on and learning gate
                    %conditions are satisfied. 
                    if curProj.learnOn
                        if isempty(curProj.learnGate)
                            curProj.update(sourceAct,curLayer.get_act()); %while on each projection, update it's weights 
                        else
%                             curProj.update(sourceAct,curLayer.get_act(),{obj.so(curProj.learnGate{1}{1}),curProj.learnGate{1}{2}});
                             if ischar(curProj.learnGate{1}{1}) %if the learnGate object is just the name, replace it with the object (backward compatibility issue)
                                curProj.learnGate{1}{1} = obj.so(curProj.learnGate{1}{1});        
                             end
                             
                            %error('abc') 
                            curProj.update(sourceAct,curLayer.get_act(),{curProj.learnGate{1}{1},curProj.learnGate{1}{2}});
%                             curProj.update(weightedAct,curLayer.get_act(),{curProj.learnGate{1}{1},curProj.learnGate{1}{2}});
                        end                        
                    end

                    
                   %%% Debug specific projections %%%
%                    if strcmp(curProj.name, 'p_conjDes_conjOut') %&& atTime(12)
%                        error()
%                    end
                   
                end
           
                
                nodes = curLayer.get_nodes();

                for node_i = 1:numel(nodes)
                    curNode = nodes{node_i};
                    nodeVals = curNode.get_vals();

                    if length(nodeVals)>1
                        weightedAct = nodeVals;
                    else
                        weightedAct = ones(1,curLayer.numUnits)*nodeVals;
                    end

                    sumExcit = sumExcit+max(weightedAct,0);
                    sumInhib = sumInhib+min(weightedAct,0);
                end

                %inhibitor = {nodeName, [start, end, dtSign], strength}; dtSign = -1,0,or 1 for neg, both, pos
                inhibitors = curLayer.get_inhibitors(); 

                
                sumOscInhib = zeros(1,curLayer.numUnits);
                for inhib_i = 1:length(inhibitors)
                    
                    
                    inhibitor = inhibitors{inhib_i};
                    if ischar(inhibitor{1}) %if the inhibitor object is just the name, replace it with the object (backward compatibility issue)
                        curLayer.inhibitors{inhib_i}{1} = obj.so(inhibitor{1});
                        inhibitors = curLayer.get_inhibitors(); 
                        inhibitor = inhibitors{inhib_i};
                    end
                    [inhibVals, inhibIndex] = max(inhibitor{1}.get_vals());
                    inhibNodeDt = inhibitor{1}.get_dt();
                    inhibNodeDt = inhibNodeDt(inhibIndex);
                    
                    intervalBounds = [inhibitor{2}(1), inhibitor{2}(2)];
                    v1 = min(intervalBounds); v2 = max(intervalBounds);
                    
%                     if strcmp(curLayer.name, 'next') && any(obj.so('n_bg').vals>0) && all(obj.so('n_bg').vals<1)                        
%                         error('abc')
%                     end

                    if ~(any(inhibVals >= v1) && any(inhibVals <= v2)) %is oscillation in the appropriate phase?
                        continue
                    elseif inhibNodeDt==0 && inhibVals==0
                        continue
                    elseif (inhibNodeDt*inhibitor{2}(3))<0
                        continue                 
                    end               
                    sumOscInhib = sumOscInhib+inhibitor{3};
                    
                    
                    
                end                              
                
                sumInhib = sumInhib+sumOscInhib;

                curLayer.update(sumExcit,sumInhib);

                if any(curLayer.act>5)
                    error('Numerical instability')
                end
               
            end
        

            for curNode = obj.nodeList

                if 0 %curNode{1}.specialUpdate>0
                    curNode{1}.specUpdate();
                else                    
                    curNode{1}.update();
                end
            end


            %flip (make temp values the new values)
            for layer_i = 1:numel(obj.layerList)
                obj.layerList{layer_i}.flip();
            end
            for proj_i = 1:numel(obj.projectionList)
                obj.projectionList{proj_i}.flip();
            end
            
            
            
            if ~isempty(obj.displayList) && mod(t,1)<(dt/2)
                %disp('DISPLAY ON')
                obj.display(); %currently broken
                
            end
                
        end
        
        function obj = disconnect(obj, sourceName, targetName)
%            error('abc')
           target = obj.so(targetName);
           
           if strcmp(targetName(1:2), 'n_')
               inputIndices = 1:numel(target.inputs);
               for input_i = 1:numel(target.inputs)  
                   input = target.inputs{input_i};
                   if strcmp(input.name, sourceName)
                       inputIndices = inputIndices(inputIndices ~= input_i);
                       target.inputs = target.inputs(inputIndices);
                       break
                   end
               end
           else
               projIndices = 1:numel(target.projections);
               for proj_i = 1:numel(target.projections)
                   proj = target.projections{proj_i};                
                   if strcmp(proj.name, ['p_' sourceName '_' targetName])                  
                       projIndices = projIndices(projIndices~=proj_i);
                       target.projections = target.projections(projIndices);
                       break
                   end
               end
           end
            
            
        end
        
        
        function obj = connect(obj,sourceName,targetName,weightConfig,norm,scaleFactor,learnType,learnRate,gates)
            numReqArgs = 8;
            
            source = obj.so(sourceName);
            target = obj.so(targetName);
            
            sourceSize = numel(source.get_act());

            targetSize = numel(target.get_act());

            newProjName = ['p_' sourceName '_' targetName];
            
            if nargin>numReqArgs
                newProj = projection(newProjName,source,target,sourceSize,targetSize,weightConfig,norm,scaleFactor,learnType,learnRate);
                for curGate = gates
                    newProj.gates = horzcat(newProj.gates,curGate);
                end
            else
                newProj = projection(newProjName,source,target,sourceSize,targetSize,weightConfig,norm,scaleFactor,learnType,learnRate);
            end
     
     
            obj.projectionList{end+1} = newProj;
            obj.accessDict(newProjName) = newProj;
            
            if strcmp(targetName(1:2),'p_')
                target.addMod(newProj);
            else
                target.addInput(newProj);
            end
        end
        
        function obj = connectNode(obj,nodeName,targetName) %nodes automatically and exclusively connect 1-to-1
            node = obj.so(nodeName);
            target = obj.so(targetName); 
            
            if node == 0
                error('No node with that name')
            end
            
            target.addNode(node);
        end
        
        function obj = log(obj,objectNames)

            for object_i = 1:numel(objectNames)                  
                obj.so(objectNames{object_i}).startLogging();
            end
           
        end
        
        function obj = logAll(obj)
            for obj_i = 1:numel(obj.layerList)
                obj.layerList{obj_i}.startLogging();
            end
            for obj_i = 1:numel(obj.projectionList)
                obj.projectionList{obj_i}.startLogging();
            end
            for obj_i = 1:numel(obj.nodeList)
                obj.nodeList{obj_i}.startLogging();
            end
        end
                
        
        function nameCell = names(obj, silent)
            if nargin < 2
                silent = false;
            end
            
            if ~silent
                cellfun(@(x) disp(x.name), obj.layerList)
                cellfun(@(x) disp(x.name), obj.nodeList)
                cellfun(@(x) disp(x.name), obj.projectionList)
            end
            nameCell = cellfun(@(x) x.name, horzcat(obj.layerList, obj.nodeList, obj.projectionList), 'UniformOutput', false);            
        end
        
        
        function obj = plot_old(obj, itemNames, nodeNames, times)
            global endTime dt
            lineTypes = {':','-','--','-.'};
            colors = {'k','b',[0,.8,.2],'r',[.26,.8,.8],'m',[.9,.9,.4], [.6,0,1],[1,.6,0]};
            numColors = 9;
            extras = 1; 
            numPlots = numel(itemNames)+numel(nodeNames)+extras;
            other_i = 0;
            if numPlots>1
                numCols = 2;
            else
                numCols = 1;
            end
            numRows = ceil(numPlots/numCols);
            scrsz = get(groot,'ScreenSize');
            figure('Position',[1 scrsz(4)/1.2 scrsz(3)/1.5 scrsz(4)/1.5]);
            
            for item_i = 1:numel(itemNames)
                itemName = itemNames{item_i};
                subplot(numRows,numCols,item_i);
                hold on
                title(itemName)
               
                if nargin == 2 %times not specified
                    switch class(obj.so(itemName))
                        case 'layer'
                            plotData = obj.so(itemName).get_actLog();
                           
                            for trace_i = 1:size(plotData,1)
                                
                                plot(max(endTime-50,dt):dt:endTime,squeeze(plotData(trace_i,:)),'LineWidth',2,...
                                    'LineStyle', lineTypes{ceil(trace_i/numColors)});
                                trace_i
                            end
                        
                        case 'projection'
                            plotData = curProj.get_weightLog('2d');
                            plot(max(endTime-50,dt):dt:endTime,plotData,'LineWidth',2)
                    end
                elseif nargin>2
               
                    switch class(obj.so(itemName))
                        case 'layer'
                            plotData = obj.so(itemName).get_actLog(times);
                            %plot(times(1):dt:times(2),plotData,'LineWidth',2)
                            for trace_i = 1:size(plotData,2)
                                
                                plot(times(1):dt:times(2),squeeze(plotData(:,trace_i)),'LineWidth',2,...
                                    'Color', colors{mod(trace_i,length(colors))+1}, 'LineStyle', lineTypes{mod(ceil(trace_i/numColors),length(lineTypes))+1});
                                
                            end
                            
                            if strcmp(itemName,'state')
                                legend1 = legend('show');
                                set(legend1,...
                                    'Position',[0.911830357142854 0.260663507109004 0.0808035714285714 0.673380726698262]);
                            end
                            %
                        case 'projection'
                            plotData = curProj.get_weightLog('2d',times);
                            plot(times(1):dt:times(2),plotData,'LineWidth',2)
                            
                    end
                    axis([times(1) times(2) min(plotData(:))-.1 max(plotData(:))+.21])
                end
                
             
            end
            
            
           
            for node_i = 1:numel(nodeNames)               
                %figure()
                subplot(numRows,numCols,item_i+other_i+node_i)
                hold on
                
                title(strrep(nodeNames{node_i},'_','\_'))
                
                %plotData = obj.so(nodeNames{node_i}).get_valsLog();
                
                if nargin == 2
                    plotData = obj.so(nodeNames{node_i}).get_valsLog();
                    plot(max(endTime-50,dt):dt:endTime,plotData,'LineWidth',2)
                    
   
                elseif nargin>2
                    plotData = obj.so(nodeNames{node_i}).get_valsLog(times);
                    plot(times(1):dt:times(2),plotData,'LineWidth',2)
                    %dbstop
                end

                axis([times(1) times(2) min(plotData(:))-.1 max(plotData(:))+.21])
                
                %plot(plotData,'LineWidth',2)                
%                 legend('show')
                
            end
            temp = 99;
        end
        
        function obj = plot(obj, itemNames, interval, special)
            %itemNames is a cell array of layers, nodes, and projections
            %whose activities, values, and weights will be plotted,
            %respectively. Interval specifies the times which will be
            %plotted either with a two-number true interval or a scalar
            %indicating how far from the end the end of the simulation to
            %start plotting. Special is an optional string that will create
            %a plot with a particular format. 
            global endTime dt t
            
            if nargin<4
                special='';
            end
            
            
            lineTypes = {':','-','--','-.'};
            colors = {'k','b',[0,.8,.2],'r',[.26,.8,.8],'m',[.9,.9,.4], [.6,0,1],[1,.6,0]};
            numColors = 9;
            numPlots = numel(itemNames);
            other_i = 0;
            if numPlots>1
                numCols = 2;
            else
                numCols = 1;
            end
            if strcmp(special,'proj')
            %make a plot showing the activity of the source and target of a
             %projection, any gate values, and overlay gate open and close
             %times
                numCols=1;
                projName=itemNames{1};
                proj=obj.so(itemNames{1}); %expecting a projection as first item name
                newList={proj.source.name, proj.target.name, proj.gates{1}{1}};
                newList=horzcat(newList,itemNames{2:end});
                itemNames=newList; 
                numPlots = numel(itemNames);
            end
            
            numRows = ceil(numPlots/numCols);
            scrsz = get(groot,'ScreenSize');
            figure('Position',[1 scrsz(4)/1.2 scrsz(3)/1.5 scrsz(4)/1.5]);
            
            if nargin==2 %if no times specified
                times=[dt,t]; %plot everything
            elseif nargin>2
                if length(interval)==2 %if two times specified
                    times=[interval(1),min(interval(2),t)]; %plot between those times
                elseif length(interval)==1 %if one number specified
                    times=[max(t-interval,dt),t]; %plot the end minus that number, to the end
                end
            end
                
            
            
            for item_i = 1:numel(itemNames)
                itemName = itemNames{item_i};
                subplot(numRows,numCols,item_i);
                hold on
                
                %font size not changing for some reason
                title(strrep(itemName,'_','\_'), 'FontSize', 18)
                                
                curItem=obj.so(itemName);
                
                if strcmp(special,'proj')

                    obj.so(projName).plotGates();
                    
                end
               
                switch class(curItem)
                    case 'layer'
                        plotData = curItem.get_actLog(times);
                        %plot(times(1):dt:times(2),plotData,'LineWidth',2)
                        for trace_i = 1:size(plotData,2)

                            plot(times(1):dt:times(2),squeeze(plotData(:,trace_i)),'LineWidth',2,...
                                'Color', colors{mod(trace_i,length(colors))+1}, 'LineStyle', lineTypes{mod(ceil(trace_i/numColors),length(lineTypes))+1});
                        end

                        if strcmp(itemName,'state')
                            legend1 = legend('show');
                            set(legend1,...
                                'Position',[0.911830357142854 0.260663507109004 0.0808035714285714 0.673380726698262]);
                        end
                        %
                    case 'projection'
                        plotData = curItem.get_weightLog('2d',times);
                        plot(times(1):dt:times(2),plotData,'LineWidth',2)
                    
                    case 'node'
                        plotData = curItem.get_valsLog(times);
                        plot(times(1):dt:times(2),plotData,'LineWidth',2)
                    
                    case 'oscNode'
                        plotData = curItem.get_valsLog(times);
                        plot(times(1):dt:times(2),plotData,'LineWidth',2)
                        
                end
                
                if strcmp(special, 'proj')
                    
                end
                
                axis([times(1) times(2) min(plotData(:))-.1 max(plotData(:))+.21])    
            end          
        end
        
        function obj = weightImage(obj, projectionName)
            weights = obj.so(projectionName).weights;
            sourceName = regexp(projectionName, '_(\w)+_', 'tokens','once');
            targetName = regexp(projectionName, '_([^_\W])+$', 'tokens', 'once');
            figure()
            imagesc(weights);
            hold on
            colorbar;
            xlabel(targetName);
            ylabel(sourceName);
            title([strrep(projectionName,'_','\_') ' weights']);
        end
            
        
        %used for debugging
        function result = conjCheck(obj,time,spec)
            global dt
            numStates = obj.so('state').numUnits;
            if nargin<3
                [~,maxUnit] = max(obj.so('conj').actLog(round(time/dt),:));
            else
                name = ['conj' spec];
                [~,maxUnit] = max(obj.so(name).actLog(round(time/dt),:));
            end
            from = ceil(maxUnit/numStates);            
            to = mod(maxUnit,numStates);
            if to == 0
                to = numStates;
            end            
            conjUnit = [from,to]
            
            switch spec
                case 'Des'
                    if obj.so('n_qOn').vals
                        [~,state] = max(obj.so('stateSim').actLog(round(time/dt),:));
                    else
                        [~,state] = max(obj.so('state').actLog(round(time/dt),:));
                    end
%                     [~,state] = max(obj.so('state').actLog(round(time/dt),:));
                    [~,next] = max(obj.so('next').actLog(round(time/dt),:));
                    correctUnit = [state,next]
                    
                case 'Obs'
                    [~,state] = max(obj.so('state').actLog(round(time/dt),:));
                    [~,last] = max(obj.so('prevState2').actLog(round(time/dt),:));
                    correctUnit = [last,state]
            end
            
        end
        
        %used for debugging
        function result = qCheck(obj,time)
            thresh = .15;
            q = obj.so('qOut').actLog(time/dt,:);
            numTrans = sum(q>thresh);
            result = zeros(numTrans,2);
            
            [~,order] = sort(q(q>thresh));
            order = fliplr(order);
            
            trans = find(q>thresh);
            trans = trans(order);
            
            
            for tr = trans
                from = ceil(tr/6);
                to = mod(tr,6);
                if to == 0
                    to = 6;
                end 
                result(tr,:) = [from,to];
            end
                
        end
         
        %used for debugging
         function result = conjcfCheck(obj,time)
            global dt
            [~,maxUnit] = max(obj.so('conjCF2').actLog(time/dt,:));
            from = ceil(maxUnit/6);            
            to = mod(maxUnit,6);
            if to == 0
                to = 6;
            end            
            result = [from,to];

         end
         
         %used for debugging
         function result = proxCheck(obj,time)
             global dt
             activations = obj.so('proxPast').actLog(time/dt,:);
             [s_act,s_i] = sort(activations);
             order = fliplr(s_i)
             activities = activations(order)
         end
             
        
         function obj = clearLogs(obj)
             global endTime dt
                for layer = obj.layerList
                    layer = layer{1};
                    if ~isempty(layer.actLog)                      
                        layer.actLog = zeros(endTime/dt,layer.numUnits);
                    end
                end
                for proj = obj.projectionList
                    proj = proj{1};
                    if ~isempty(proj.weightLog)
                        proj.weightLog = zeros(size(proj.weights,1),size(proj.weights,2),endTime/dt);
                    end
                    if ~isempty(proj.traceLog)
                        proj.traceLog = zeros(endTime/dt,size(proj.weights,2));
                    end
                end
                for node = obj.nodeList
                    node = node{1};
                    if ~isempty(node.valsLog)
                        node.valsLog = zeros(endTime/dt,numel(node.vals));
                    end
                end
         end
         
         function obj = reset(obj)
             obj.clearLogs();
             for layer = obj.layerList
                 layer = layer{1};
                 layer.act = layer.act*0;
                 layer.actTemp = layer.actTemp*0;
                 %traces? currently no layer traces used
             end
             
             for proj = obj.projectionList
                 proj = proj{1};
                 if ~isempty(proj.timeTrace)
                    proj.timeTrace = proj.timeTrace*0;
                 end
             end
         end
         
         function obj = set_display(obj,actPlotList, envPlotList)
             for itemName = actPlotList
                 obj.displayList{1}{end+1} = obj.so(itemName{1});
             end
             for itemName = envPlotList
                 obj.displayList{2}{end+1} = obj.so(itemName{1});
             end
         end
         
         
         function obj = exportPlotData(obj,itemList,times, filename, savedir)
             global dt endTime
             oldDir = pwd;
             newDir = ['/data/hammer/space0/model/nz_goals/' savedir];
             if ~exist(newDir, 'dir')
                 mkdir(newDir)
             end
             cd(newDir);
             
             if strcmp(itemList, 'all')
                 itemList = obj.names(true);
             end
             
             if isempty(times)
                 times=[dt endTime];
             end
             
             %write text file
             fileID = fopen('saveDetails.txt', 'w');
             fprintf(fileID, 'Date ran: %s\n', date);
             fprintf(fileID, 'Times: %s\n', num2str(times));
             fprintf(fileID, 'Items:');
             fprintf(fileID, '%s, ', itemList{:});
             fclose(fileID);
             
             
%              error('abc');
             
             for itemName = itemList
                 itemName = itemName{1};
                 fullFileName = [itemName '_' filename '.csv'];
                 
                 if strcmp(itemName(1:2), 'p_')
                     fileData = obj.so(itemName).weightLog;
                     allTimes = times(1)/dt:times(2)/dt;
                     nrow = size(fileData, 1); ncol = size(fileData, 2);
                     fileData = fileData(:,:,allTimes);
                     fileData = permute(fileData, [3, 1, 2]);
                     fileData = reshape(fileData, length(allTimes), nrow*ncol);
                     
                     error('abc')
                 
                     
                 elseif strcmp(itemName(1:2), 'n_')
                     
                    fileData = obj.so(itemName).valsLog;                 
                    fileData = fileData(times(1)/dt:times(2)/dt,:);
                     
                 else
                     
                    fileData = obj.so(itemName).actLog;                 
                    fileData = fileData(times(1)/dt:times(2)/dt,:);
                 end
               
                 
                 csvwrite(fullFileName,fileData);
                 cd(oldDir);
             end
         end
         
         
         
         
         function obj = testPlot(obj)
             global t dt
             data=ones(1,t/dt) * .5;
             plot(data)
         end
         
         function obj = gatePlot(obj, projName, gateName, interval)
             %make a plot showing the activity of the source and target of a
             %projection, any gate values, and overlay gate open and close
             %times
            global endTime dt t

            if nargin == 3
                interval=[dt, t];
            end
            
            lineTypes = {':','-','--','-.'};
            colors = {'k','b',[0,.8,.2],'r',[.26,.8,.8],'m',[.9,.9,.4], [.6,0,1],[1,.6,0]};
            numColors = numel(colors);
            
            
            numCols=1;
            proj=obj.so(projName); %expecting a projection as first item name
            itemNames={proj.source.name, proj.target.name, gateName};
            if strcmp(itemNames{1}, itemNames{2})
                itemNames(2)=[];
            end
            
            if ~isempty(proj.learnGate)              
                if strcmp(proj.learnGate{1}{1}.name, gateName)
                    plotTraces = true;
                else
                    plotTraces = false;
                end
            else
                plotTraces = false;
            end
            
            numPlots = numel(itemNames);            
            numRows = ceil(numPlots/numCols);
            scrsz = get(groot,'ScreenSize');
            figure('Position',[1 scrsz(4)/1.2 scrsz(3)/1.5 scrsz(4)/1.5]);
            
            if nargin==2 %if no times specified
                times=[dt,t]; %plot everything
            elseif nargin>2
                if length(interval)==2 %if two times specified
                    times=[interval(1),min(interval(2),t)]; %plot between those times
                elseif length(interval)==1 %if one number specified
                    times=[max(t-interval,dt),t]; %plot the end minus that number, to the end
                end
            end
            
            
            %get gate open times
            [openTimes, closeTimes] = getGateIntervals(obj, projName, gateName, interval);
            hold on
            openx=zeros(2,numel(openTimes));
            openy=zeros(2,numel(openTimes));
            openy(2,:)=1;
            for time_i = 1:numel(openTimes)
                openx(:,time_i) = [openTimes(time_i),openTimes(time_i),]';
            end
             
            closex=zeros(2,numel(closeTimes));
            closey=zeros(2,numel(closeTimes));
            closey(2,:)=1;
            for time_i = 1:numel(closeTimes)
                closex(:,time_i) = [closeTimes(time_i),closeTimes(time_i),]';
            end
            
            %get trace info
            if plotTraces
                traceData = proj.get_traceLog(times);
            end
            
            for item_i = 1:numel(itemNames)
                itemName = itemNames{item_i};
                subplot(numRows,numCols,item_i);
                hold on
                title(strrep(itemName,'_','\_'))
                
                curItem=obj.so(itemName);
   
               
                switch class(curItem)
                    case 'layer'
                        plotData = curItem.get_actLog(times);
                        %plot(times(1):dt:times(2),plotData,'LineWidth',2)
                        hold on
                        for line_i = 1:size(plotData,2)

                            plot(times(1):dt:times(2),squeeze(plotData(:,line_i)),'LineWidth',2,...
                                'Color', colors{mod(line_i,length(colors))+1}, 'LineStyle', lineTypes{mod(ceil(line_i/numColors),length(lineTypes))+1});
                            
                            if plotTraces
                                 plot(times(1):dt:times(2),squeeze(traceData(:,line_i)),'LineWidth',2,...
                                'Color', colors{mod(line_i,length(colors))+1}, 'LineStyle', ':');
                            end
                            
                        end                 
                    case 'oscNode'
                        hold on
                        plotData = curItem.get_valsLog(times);
                        plot(times(1):dt:times(2),plotData,'LineWidth',2)
                        
                        openy(1,:)=-1;
                        closey(1,:)=-1;             
                end
                
                %add gate lines
                line(openx,openy, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
                line(closex,closey, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
                
                axis([times(1) times(2) min(plotData(:))-.1 max(plotData(:))+.21])    
            end      
%             error('abc')
         end
         
         
        
        
         
         function [openTimes, closeTimes] = getGateIntervals(obj, projName, gateName, interval)
            global t dt
           
            if length(interval)==2 %if two times specified
                times=[interval(1),min(interval(2),t)]; %plot between those times
            elseif length(interval)==1 %if one number specified
                times=[max(t-interval,dt),t]; %plot the end minus that number, to the end
            end
             
            proj = obj.so(projName);
            possGates = horzcat(proj.gates, proj.learnGate);
            for gate_i = 1:numel(possGates)
                tryGate = possGates{gate_i};
                if strcmp(tryGate{1}.name, gateName)
                    gate=tryGate{1};
                    gateConditions = tryGate{2};
                    if length(gateConditions)<3
                        gateConditions(3) = 0;
                    end
                    break
                end
            end
            gateVals=gate.get_valsLog(times)';
            gateDvdt=gate.get_dvdtLog(times)';      

            v1 = min(gateConditions(1:2)); v2 = max(gateConditions(1:2));
            
            above=gateVals > v1;
            below=gateVals < v2;
            if gateConditions(3) ~= 0
                correctDt=(gateDvdt.*gateConditions(3)) > 0;
            else
                correctDt=ones(size(above));
            end
%             gateOpen = (above+below+correctDt)==3;
            gateOpen = (above & below & correctDt);
            
          
            indices=1:numel(gateOpen);
           
            shiftForward=horzcat(99,gateOpen(1:end-1));
            shiftBackward=horzcat(gateOpen(2:end),99);
            
            openTimes=indices(gateOpen==1 & shiftForward==0);
            closeTimes=indices(gateOpen==1 & shiftBackward==0);
            
            openTimes=openTimes.*dt + times(1);
            closeTimes=closeTimes.*dt + times(1);
           
%             error('abc')               
         end
             
       
             
                 
                 
         
         
        function obj = display(obj)
            global t
            
            %activity plots
            if length(obj.displayList{1}) > 1 
               
                numPlots = numel(obj.displayList{1});
                if numPlots>1
                    numCols = 2;
                else
                    numCols = 1;
                end
                numRows = ceil(numPlots/numCols);
                hold off
                figure(37)
                %disp('DISPLAY CALLED (line 462)')
                hold on
    %             scrsz = get(groot,'ScreenSize');
    %             figure('Position',[1 scrsz(4)/1.2 scrsz(3)/1.5 scrsz(4)/1.5]);

                for item_i = 1:numPlots
                    item = obj.displayList{1}{item_i};
                    g=subplot(numRows,numCols,item_i);

                    if strcmp(class(item),'layer')
                        imagesc(item.act,[0,1])
    %                     p=get(g,'position');
    %                     p(4)=p(4)*.1;
    %                     set(g, 'position', p);
    % %                     error()
                    elseif strcmp(class(item),'projection')
                        imagesc(item.weights,[0,1])
                    end
                    title([item.name '  t = ' num2str(t)])
                end
            end
            
            if ~isempty(obj.displayList{2})
                numPlots = numel(obj.displayList{2});
                eDims=[4,4];
                eSize=eDims(1)*eDims(2);
                itemIndicesToRemove=[];
%                 
%                 for plot_i = 1:numPlots
%                     item = obj.displayList{2}{plot_i};
%                     figure(70+plot_i)
%                     
%                     plotData = reshape(item.act, eDims(1), eDims(2));
%                     imagesc(plotData)
%                     
%                 end
                
                figure(38)
                hold on


                for item_i = 1:numPlots
                    item = obj.displayList{2}{item_i};
                    if strcmp(class(item),'layer')
                        if numel(item.act) ~= eDims(1)*eDims(2)
                            warning(['Environmental plot object ' item.name ' input size not correct'])
                            itemIndicesToRemove(end+1)=item_i;                                                     
                            continue
                        end
                        
                        g=subplot(numRows,numCols,item_i);
                        plotData = reshape(item.act, eDims(1), eDims(2));
                        imagesc(plotData)
    %                     p=get(g,'position');
    %                     p(4)=p(4)*.1;
    %                     set(g, 'position', p);
    % %                     error()
                    elseif strcmp(class(item),'projection')
                        if numel(item.weights) ~= eSize^2
                            warning(['Environmental plot object ' item.name ' input size not correct'])
                            itemIndicesToRemove(end+1)=item_i;
                            continue
                        end
                        item = obj.displayList{2}{item_i};
                        g=subplot(numRows,numCols,item_i);
                        plotData = reshape(item.weights, eSize, eSize);
                        imagesc(plotData)
    %                     p=get(g,'position');
    %                     p(4)=p(4)*.1;
    %                     set(g, 'position', p);
    % %                     error()
                    end
                    title([item.name '  t = ' num2str(t)])
                end
                
                obj.displayList{2}(itemIndicesToRemove)=[];
                
            end
        end
        
        function values = getPatterns(obj, layerName, times)
            global t dt
            %expecting time to be an nx2 matrix of time intervals
            
            if strcmp(layerName(1:2), 'n_') %if 'layerName actually refers to a node (unusual)
                node = obj.so(layerName);
                values = nan(size(times,1), length(node.vals));
                for t_i = 1:size(times,1)
                    interval = times(t_i,:)/dt;
                    logVals = node.valsLog(interval(1):interval(2),:);
                    values(t_i,:) = mean(logVals, 1);
                end
            else
            
                layer = obj.so(layerName);
                values = nan(size(times,1), layer.numUnits);
                for t_i = 1:size(times,1)
                    interval = times(t_i,:)/dt;
                    logVals = layer.actLog(interval(1):interval(2),:);
                    values(t_i,:) = mean(logVals, 1);
                end       
            end
        end
        
        function obj = exportPlotData(obj,itemList,times, filename, savedir)
             global dt endTime
             oldDir = pwd;
             newDir = ['/data/hammer/space0/model/nz_goals/' savedir];
             if ~exist(newDir, 'dir')
                 mkdir(newDir)
             end
             cd(newDir);
             
             if strcmp(itemList, 'all')
                 itemList = obj.names(true);
             end
             
             if isempty(times)
                 times=[dt endTime];
             end
             
             %write text file
             fileID = fopen('saveDetails.txt', 'w');
             fprintf(fileID, 'Date ran: %s\n', date);
             fprintf(fileID, 'Times: %s\n', num2str(times));
             fprintf(fileID, 'Items:');
             fprintf(fileID, '%s, ', itemList{:});
             fclose(fileID);
                          
             for item_i = 1:numel(itemList)
                 itemName = itemList{item_i};
                 fullFileName = [itemName filename '.csv'];
                 
                 if strcmp(itemName(1:2), 'p_')
                     
                     allTimes = times(1)/dt:times(2)/dt;
                     if obj.so(itemName).learnRate > 0
                        fileData = obj.so(itemName).weightLog;
                     else                    
                         fileData = repmat(obj.so(itemName).weights, [1,1,numel(allTimes)]);
                     end                            
                         
                     
                     nrow = size(fileData, 1); ncol = size(fileData, 2);
                     fileData = fileData(:,:,allTimes);
                     fileData = permute(fileData, [3, 1, 2]);
                     fileData = reshape(fileData, length(allTimes), nrow*ncol);
                     
             
                 elseif strcmp(itemName(1:2), 'n_')
                     
                    fileData = obj.so(itemName).valsLog;                 
                    fileData = fileData(times(1)/dt:times(2)/dt,:);
                     
                 else
                     
                    fileData = obj.so(itemName).actLog;                 
                    fileData = fileData(times(1)/dt:times(2)/dt,:);
                 end
                           
                 csvwrite(fullFileName,fileData);                 
             end
             cd(oldDir);
        end      

                         
    end   
end