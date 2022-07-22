classdef stateSpace < handle
    properties
        transitions
        curState
        stateLog = []
        adjacencyMat = [];
        numStates = 0;
        timeInState = [];
        rewards = [];
        numRewards=0
    end
    
    methods
        function obj = stateSpace(transitions, curState, rewards)
            obj.transitions = transitions;
            obj.curState = curState;
            obj.numStates = numel(unique(obj.transitions(:,1)));
            obj.timeInState = zeros(1,obj.numStates);
            
            obj.adjacencyMat = zeros(obj.numStates);
            for t_i  = 1:size(obj.transitions,1)
                obj.adjacencyMat(obj.transitions(t_i,1),obj.transitions(t_i,2)) = 1;
            end
            
            obj.stateLog(end+1,:) = [curState,0];
            
            if nargin > 2
               obj.rewards = rewards; %each row = [reward, state]
               obj.numRewards = numel(unique(obj.rewards(:,1)));
            end
            
        end
        
        function rewardVec = getReward(obj)
            rewardVec = zeros(1, obj.numRewards);
            currentState = find(obj.timeInState > 2); %only reward if been in state for two timesteps
            if isempty(currentState)
                currentState = 0;
            end
            rewardVec(obj.rewards(obj.rewards(:,2)==currentState,1))=1;
%             if atTime(20)
%                 error('abc')
%             end
           
        end
        
        
        function obj = update(obj)
            global dt
            newTime = obj.timeInState(obj.curState) + dt;
            obj.timeInState = zeros(1,obj.numStates);
            obj.timeInState(obj.curState) = newTime;
        end
        
        
        
        function value = get_loc(obj)
            value = obj.curState;
        end
        
        function value = get_locVec(obj)
            value = zeros(1,obj.numStates);
            value(obj.curState) = 1;
        end
        
        function obj = set_loc(obj,location)
            obj.curState = location;
        end
        
        
        function value = get_numStates(obj)
            value = obj.numStates;
        end
        
        function value = get_numActions(obj)
            value = max(obj.transitions(:,3));
        end
        
        function value = get_numRewards(obj)
            value = obj.numRewards;
        end
        
        function moved = rwalk(obj)
            global t
            %randomly move
            possTrans = obj.transitions(obj.transitions(:,1) == obj.curState,:);
            newLoc = possTrans(randi(size(possTrans,1)),2);
            obj.set_loc(newLoc);
            obj.get_loc();
            obj.stateLog(end+1,:) = [newLoc,t];
            moved=1;
        end
        
        %like rwalk, but return action as well
        function [moved, action] = rwalk2(obj)
            global t
            %randomly move
            possTrans = obj.transitions(obj.transitions(:,1) == obj.curState,:);
            moveResults = possTrans(randi(size(possTrans,1)),2:3);
            newLoc = moveResults(1);
            action = moveResults(2);
            obj.set_loc(newLoc);
            obj.get_loc();
            obj.stateLog(end+1,:) = [newLoc,t];
            moved=1;
            
        end
        
        function moved = nwalk(obj,nextAct) %takes in the 'next' activities and moves to the max one if possible
            global t
            actThresh=.6;
            [val,maxNext] = max(nextAct);
            if val>actThresh
                tr = obj.transitions;
                if any(tr(tr(:,1) == obj.curState & tr(:,2) == maxNext))
                    obj.set_loc(maxNext);
                end
                obj.stateLog(end+1,:) = [maxNext,t];
                moved=1;
            else
                moved=0;
            end
        end
        
        function moved = mwalk(obj,input)
            %move agent based on a motor command. The command (input) is
            %either a ector where one unit=1, in which case the index specifies
            %the action, or it is a scalar specifying the action directly
            global t
            if length(input)>1
                [val,action] = max(input);
            else
                action = input;
                val = 1;
            end
            if val>.9999
                transFromState_i = find(obj.transitions(:,1) == obj.curState);
                transWithAct_i = find(obj.transitions(:,3) == action);
                newState = obj.transitions(intersect(transFromState_i,transWithAct_i),2);
                if ~isempty(newState)
                    obj.set_loc(newState);
                    moved = 1;
                    obj.stateLog(end+1,:) = [newState,t];
                else
%                     warning('Cannot take that action here!')
                    moved = -1;
                end
                
            else
                moved = 0;
            end
        end
        
    end
end
        