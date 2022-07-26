classdef layer < handle
    
    
% Layer object consisting of model neurons
%     
% 
% Main Properties
%    act - vector of length n for n units representing unit activity
%    actTemp - temporary representation of updated activity
%    noiseGain - controls the level of noise in the neural activity
%    projections - cell array of projections projecting onto this layer
%    nodes - cell array of nodes (non-neural input) projecting onto this layer
%    inhibitors - cell array of nodes or layers globally inhibiting this
%       layer
%    actFunc - name of activation function
%    timeConst - time constant used in updating activity via actFunc
%    decayRate - decay rate used in updating activity via actFunc
%    logOn - controls whether log is being updated
%    actLog - log of activity
%    dadt - vector of activity time derivative


% 
% Main Methods
%     layer - sets up layer with specified parameters
%      update - updates layer activity using the total excitatory and
%      inhibitory activity coming into the layer (calculated in the network
%      object update method).
    
    properties
        name
        numUnits
        act
        actTemp
        noiseGain=0;
        bias=0; %not currently used
        
        projections={} 
        nodes={}
        inhibitors={}
        
        actFunc=''
        timeConst
        decayRate
          
        logOn=false
        actLog=[]        
        dadt=[];
        
        %obsolete properties
        %baseline=0;
        %artTrace=[];
      
        %traceType='' %does it ever make sense to have traces in the layer instead of the projection?
        %timeTrace=[]
        %traceLog=[]
        %gates={}
        
    end
    
    methods
        function obj=layer(name,numUnits,actFunc,timeConst,decayRate,traceType,inhibitors)
            obj.name=name;
            obj.numUnits=numUnits;
            obj.actFunc=actFunc;
            obj.timeConst=timeConst;
            obj.decayRate=decayRate;
            
            obj.projections={};
            obj.act=zeros(1,numUnits);
            obj.actTemp=zeros(1,numUnits);
            obj.dadt=obj.act;
            
            
            %obsolete
%             obj.baseline=0;            
            
            if nargin>5 %was 6 when layer-level inhibitors were still possible
                obj.inhibitors=inhibitors; %{{nodename,inhibitionTiming,inhibitionStrength}};
            end
        end
        

        
        function obj=update(obj,totalExcit,totalInhib)
            global dt t
      
            switch obj.actFunc
                
                case 'shunt'

                      obj.actTemp = obj.act + dt*((-obj.decayRate*obj.act + (1-obj.act).*totalExcit...
                        + obj.act.*totalInhib + randn(1,obj.numUnits)*obj.noiseGain+obj.bias)/obj.timeConst);

               
                case 'latInhib' %linear, not shunting, inhibition (what is first min doing?)

%                     obj.actTemp = max(obj.act + dt*(( min(-obj.decayRate*obj.act,0) + (1-obj.act).*totalExcit...
%                         + totalInhib + randn(1,obj.numUnits)*obj.noiseGain+obj.bias)/obj.timeConst),0);

                    obj.actTemp = max(obj.act + dt*(( -obj.decayRate*obj.act + (1-obj.act).*totalExcit...
                        + totalInhib + randn(1,obj.numUnits)*obj.noiseGain+obj.bias)/obj.timeConst),0);
                    

                case 'proxMap' %makes units easier to inhibit by making the inhibition factor=1 if act>.2, ratehr than the factor=act
                    obj.actTemp = max(obj.act + dt*(( min(-obj.decayRate*obj.act,0) + (1-obj.act).*totalExcit...
                        + (min(obj.act,.2)*5).*totalInhib + randn(1,obj.numUnits)*obj.noiseGain+obj.bias)/obj.timeConst),0);
                                                                 
                    
                case 'queue' % a unit only decays when its activity is <=.15
                    decayThresh=.15;
                    obj.actTemp = obj.act + dt*((-obj.decayRate*(obj.act<=decayThresh).*(obj.act>0) + (1-obj.act).*totalExcit...
                        + obj.act.*totalInhib + randn(1,obj.numUnits)*obj.noiseGain+obj.bias)/obj.timeConst);
                    
                case 'dummy'
                    obj.actTemp = totalExcit + totalInhib;
                    
                otherwise
                    error('Invalid activation function')
            end
            

            
            if obj.logOn
                obj.actLog(round(t/dt),:)=obj.act;
                               
            end
            
            obj.actTemp=max(obj.actTemp,0);            
            obj.dadt=(obj.actTemp-obj.act)/dt;                   
        end
        
        function obj=addInput(obj,proj)
            obj.projections{end+1}=proj;
        end
        
        function obj=flip(obj)
            obj.act=obj.actTemp;
        end
        
        function value=get_projections(obj)
            value=obj.projections;
        end
        
        function value=get_nodes(obj)
            value=obj.nodes;
        end
        
        function obj=set_inhibitors(obj,inhib_array)
            
            for inhibitor=inhib_array
%                 if strcmp(obj.name,'adjFuture')
%                     error('debug')
%                 end
                inhibitor=inhibitor{1};
                if length(inhibitor)<3 %if missing inhibition strength
                    inhibitor{3}=-20; %make it -10
                end
                obj.inhibitors{end+1}=inhibitor;               
            end
        end
        
        function value=get_inhibitors(obj)
            value=obj.inhibitors;
        end
        
        function obj=addNode(obj,node)
            obj.nodes{end+1}=node;
        end
        
        function i=maxUnit(obj, time)
            global dt
            if nargin<2
                [~,i]=max(obj.act);
            else
                [~,i]=max(obj.actLog(round(time/dt),:));
            end
        end
        
        
        function value=get_act(obj)
            value=obj.act;
        end
        
        function obj=startLogging(obj)
            global endTime dt
            obj.logOn=true;
            obj.actLog=zeros(endTime/dt,obj.numUnits);

        end
        
        function value=get_actLog(obj,times) %times=[startTime,stopTime]
            global dt
            if nargin>1
                value=obj.actLog(round(times(1)/dt):round(times(2)/dt),:);
            else
                value=obj.actLog;
            end
        end
        
        function obj=clearAct(obj)
            obj.act=obj.act*0;
            obj.actTemp=obj.act;
        end           
         
    end       
            
end
