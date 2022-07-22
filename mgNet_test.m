%sets up network for small 6-state environment simulation
% Testing modifications to 'find next state' components after April '17 cleanup.
% among other things, changed proxPast learning rule to be oscillatory and
% cleaned up plot method. 

NEW_NETWORK=1;
if NEW_NETWORK
    clear all
    NEW_NETWORK=1;
else
    n.reset()
end

CLOSE_ALL=1; %closes all figures right before creating new ones (after the run)
SET_WEIGHTS=0; %initialize weights to correct values
rootdir='/data/hammer/space0/model/nz_goals/';
addpath(rootdir);
cd(rootdir);

global dt endTime t 
dt=.05;

START_STATE = 1;

%set up network and environment
if NEW_NETWORK
    
    
    load('/data/hammer/space0/model/nz_goals/stateSpaces/mgEnv.mat'); %script that sets up variable
    load('/data/hammer/space0/model/nz_goals/savedWeights/mgConjWeights.mat');
    e=stateSpace(mgEnv,START_STATE); %(transitionMatrix, starting location)
    
    numStates =  e.get_numStates();
    numActions = e.get_numActions();
    numInps=numStates;
 
    n=network();
    
    %LAYERS AND NODES
    %addLayer(name,{numUnits,actFunc,timeConst,decayRate,traceType(optional)})
    %layer.set_inhibitors(nodeName, [v1, v2, dt-sign]); inhibition will
    %only occur when node value is between v1 and v2 with the appropriate
    %dt sign (0 = positive or negative) 
    
    
  
    %=========================
    %== Choosing the next state ==
    %=========================
    
    n.addLayer('state',{numStates,'shunt',.5,1}); %was tc=1, dr=.5
    n.addLayer('goalOut',{numStates,'shunt',1,1});
    
    n.addLayer('proxPast',{numStates,'proxMap',1,1}); 
    
    n.so('proxPast').set_inhibitors({{'n_novInhib',[.8,1,-1],-8},{'n_qOsc',[-.6,-.3,1]},{'n_subOsc',[-.2,.2,0]}});
    
  
    n.addLayer('adjFuture',{numStates,'latInhib',.2,1}); 
    n.so('adjFuture').set_inhibitors({{'n_novInhib',[.7,1,-1]},{'n_qOsc',[-.9,-.3,1]}});
    
    n.addLayer('next',{numStates,'latInhib',1,.1});
    n.so('next').noiseGain=.025;
%     n.so('next').set_inhibitors({{'n_novInhib',[.15,.05,-1]}, {'n_qOsc',[-.4,.1,1]}});   %{'n_qOsc',[.01, .4,1]}});
    n.so('next').set_inhibitors({{'n_novInhib',[.5,.05,-1]}, {'n_qOsc',[-.4,.1,1]}});   %{'n_qOsc',[.01, .4,1]}});

    
    %=========================
    %== Monitoring the past ==
    %=========================
    
    n.addLayer('prevState1',{numStates,'shunt',4,.5}); %tc=10
    n.addLayer('prevState2',{numStates,'latInhib',1,.001});
    n.so('prevState2').set_inhibitors({{'n_novInhib',[.7,1,0]}});
     
    n.addLayer('prevMotor1',{numActions,'shunt',4,.5}); %tc=10
    n.addLayer('prevMotor2',{numActions,'latInhib',1,.001});
    n.so('prevMotor2').set_inhibitors({{'n_novInhib',[.6,1,0]}});
    
    n.addLayer('conjObs', {numStates^2, 'latInhib',1.5,1});%tc=.5
    n.so('conjObs').set_inhibitors({{'n_novInhib',[.5,.8,0]}});
    
    %=========================
    %== Determining Action ==
    %=========================
    n.addLayer('conjDes', {numStates^2, 'latInhib',1.5,1});
    n.so('conjDes').set_inhibitors({{'n_novInhib',[.1,1,-1]},{'n_qOsc',[.1, .8, 1]}});
    
    n.addLayer('conObs', {numStates^2, 'latInhib',1.5,1});
    n.so('conjObs').set_inhibitors({{'n_novInhib',[.1,1,-1]}});
    
    n.addLayer('conjOut', {numStates^2, 'latInhib',1.5,1});
    n.so('conjOut').set_inhibitors({{'n_subOsc',[-.2,.2,0]},{'n_qOsc', [-.2,.2,1]}});%{'n_subOsc',[.1, .8, 1]}
    
    n.addLayer('motorIn',{numActions,'latInhib',1,1}); %6 actions: left->mid, left->right, mid->left, etc.; TC, decay was .2,2 9/29
    n.so('motorIn').set_inhibitors({{'n_subOsc',[-.5,.5,0]}});
    
    n.addLayer('motorOut',{numActions,'latInhib',.5,.2});
    
    %============================
    %== Storing a Transition Plan ==
    %=============================
    n.addLayer('stateSim',{numStates,'latInhib',.7,.01});
    n.so('stateSim').set_inhibitors({{'n_novInhib',[.6,1,-1]},{'n_qOsc',[-.7,-1,-1]}});

    n.addLayer('qIn',{numStates^2, 'latInhib', .2,1});
    n.so('qIn').set_inhibitors({{'n_qOsc',[.05, .15, 1]}}); %[-.6, -.7, -1]}});
    
    n.addLayer('qStore',{numStates^2, 'queue', 1,1});
    
    n.addLayer('qOut', {numStates^2, 'latInhib', .2,0});
    n.so('qOut').set_inhibitors({{'n_qOsc', [-.2,.2,1]}});
    
      %============================
    %== Nodes ==
    %=============================
    
    %External Inputs
    n.addNode('n_env',numStates);
    n.setNode('n_env',e.get_locVec);
    n.connectNode('n_env','state');

    n.addNode('n_goalOut',numStates);
    
    %Action
    bgFunc=@(vals,in) (any(vals>0)*max((vals-(dt/5)),0))+(all(vals<=0).*(in>=.6).*(in==max(in)));
    n.addNode('n_bg',numActions);
    n.so('n_bg').updateFunc=bgFunc;
    n.so('n_bg').inputs={n.so('motorOut')};
    
    
    n.addNode('n_novInhib',1);
    n.so('n_novInhib').updateFunc = @(x) max(x - (dt/2), 0);
    
    %Oscillations
    n.addOscNode('n_subOsc',2);
    n.addOscNode('n_qOsc',1); %1.5

    %Control signals  
    n.addNode('n_qOn',1); % 0=no cue, 1=storing plan, 2=reading plan
    
    %PROJECTIONS
    %(obj,sourceName,targetName,weightConfig,norm,weightFactor,learnType,learnRate,gates)
    %gate format = {'name', [val1, val2]} gate open between val1 & val2
    
    
    %=========================
    %== Choosing the next state ==
    %=========================
    
    n.connect('state','adjFuture','one-to-one',0,2,'adj_future',1,{{'n_qOn',[0 0]}});
    n.connect('state','proxPast','one-to-one',0,4,'null',1,{{'n_subOsc',[.2,1]}});
    
    n.connectNode('n_goalOut','goalOut'); 
    n.connect('goalOut','proxPast','one-to-one',0,2,'null',0, {{'n_subOsc',[-.2,-1]}}); 
   
    n.connect('proxPast','proxPast','one-to-one',2,0,'prox_past',1); 
    n.so('p_proxPast_proxPast').learnGate={{'n_subOsc',[.2,1]}};
    
    n.connect('proxPast', 'next', 'one-to-one', 0,2,'null',0, {{'n_subOsc',[-1,-.2,1]}});
    
    
    n.connect('next','next','impose',0,2,'null',0);
    
    n.so('p_proxPast_next').channels={n.so('adjFuture'),.3};
 
    
     %=========================
    %== Monitoring the past ==
    %=========================
    n.connect('state','prevState1','one-to-one',0,0,'null',0);
    n.connect('prevState1','prevState2','one-to-one',0,1,'null',0);
    n.connect('prevState2','prevState2','impose2',0,1,'null',0); %inhibition was 2; was impost
    
    n.connectNode('n_bg','prevMotor1');
    n.connect('prevMotor1','prevMotor2','one-to-one',0,1,'null',0);%,{{'n_novInhib',[.75,1]}});%no,{{'n_novGate',-1}}); %connection made during simulation as a hack for first state   
    n.connect('prevMotor2','prevMotor2','impose2',0,2,'null',0); %i
    n.connect('prevMotor2','motorIn','impose',0,3,'null',0,{{'n_subOsc', [.2, 1]}});%,{{'n_motorGate',1}});
    
    
    n.connect('state', 'conjObs', 'null', 0, 0, 'null', 0)
    n.connect('prevState2', 'conjObs', 'null', 0, 0, 'null', 0)
    n.connect('conjObs','conjObs','impose',0,2,'null',0);
    n.connect('conjObs' ,'conjOut', 'impose2', 0, 1, 'null', 0, {{'n_subOsc', [.2, 1]}})
    %=========================
    %== Determining Action ==
    %=========================
    
     n.connect('state','conjDes','unif',0,1,'null',0,{{'n_qOn',[0,0]}});    
   

    n.connect('next','conjDes','unif',0,1,'null',0); 
    
    
    n.connect('conjDes','conjOut','impose2',0,1,'null',0, {{'n_subOsc', [-.2,-1]}, {'n_qOn', [0,0]}});
    n.connect('conjDes','conjDes','impose',0,2,'null',0);
    n.connect('conjOut','conjOut','impose',0,2,'null',0);
    n.connect('conjOut', 'motorIn','null',0,0,'oscHebb',1, {{'n_subOsc',[-.2,-1]}}); %learning rate was 20
    n.so('p_conjOut_motorIn').learnGate={{'n_subOsc',[.2,1]}};

    n.connect('motorIn','motorOut','impose',0,6,'null',0,{{'n_subOsc', [-.4,-1]}});%,{{'n_motorGate',1}});
   

     %============================
    %== Storing a Transition Plan ==
    %=============================
    n.connect('state','stateSim','one-to-one',0,1, 'null', 1, {{'n_qOn',[0 0]}}); %state projects to stateSim when queue off

    n.connect('next','stateSim', 'one-to-one',0, 6,'null',0,{{'n_qOn',[.99 1.01]},{'n_qOsc',[-.9,-.6,1]}});

    n.connect('stateSim', 'stateSim', 'impose', 0, 1, 'null', 1);
    n.connect('stateSim', 'adjFuture', 'one-to-one', 0, 2, 'adj_future', 1, {{'n_qOn', [.99 1.01]}});
   

    n.connect('stateSim', 'conjDes', 'unif', 0, 1, 'null', 0, {{'n_qOn', [.99 1.01]}});  
    n.so('p_stateSim_adjFuture').learnGate={{'n_qOn',[0 0]}};

    n.connect('conjDes', 'qIn', 'one-to-one', 0, 8, 'null', 1,{{'n_qOn', [.99 1.01]},{'n_qOsc',[-.2,0,1]}});
    n.connect('qIn', 'qStore', 'one-to-one', 0, 12, 'null', 1);
    n.connect('qStore', 'qIn', 'all-to-all', 0, -2.5, 'null', 1);
    
    n.connect('qStore','qOut','one-to-one',0,2,'null',1, {{'n_qOn', [1.01,2.1,0]}, {'n_qOsc', [.8,1,0]}});
    n.connect('qOut','qOut','impose',0,2,'null',0);
    
    n.connect('qOut','qStore','one-to-one',0,-5,'null',1,{{'n_qOsc', [-.7,-1,0]}});
    
    n.connect('qOut','conjOut','one-to-one',0,2,'null',1);%,{{'n_qOsc', [.3,-.3, -1]}});
end

% n.log({'state','adjFuture','p_state_adjFuture','n_novInhib','proxPast','p_proxPast_proxPast', 'n_subOsc', 'goalOut'});

%plotItems = {'state','stateSim','adjFuture', 'conjDes', 'proxPast','next', 'n_novInhib', 'n_qOsc','n_subOsc','goalOut','qIn','qStore','conjOut'};
%otherItems = {'p_stateSim_adjFuture', 'p_state_adjFuture', 'p_proxPast_next','p_proxPast_proxPast'};
%n.log(horzcat(plotItems, otherItems));
n.logAll();


 %============================
%== Set Weights ==
%=============================
% hardcoded weights
n.so('p_state_conjDes').set_weights(conjFromWeights);
n.so('p_stateSim_conjDes').set_weights(conjFromWeights);
n.so('p_next_conjDes').set_weights(conjToWeights);

n.so('p_prevState2_conjObs').set_weights(conjFromWeights);
n.so('p_state_conjObs').set_weights(conjToWeights);
% (typically) learned weights
if SET_WEIGHTS
    n.so('p_proxPast_proxPast').set_weights(e.adjacencyMat+eye(e.numStates));
    n.so('p_state_adjFuture').set_weights(e.adjacencyMat+eye(e.numStates).*2);
    n.so('p_stateSim_adjFuture').set_weights(e.adjacencyMat+eye(e.numStates).*2);
    
    
%     n.so('p_conjOut_motorIn').set_weights(conjActWeights);
%     n.so('p_conjOut_motorIn').learnRate=0;
    
    n.so('p_state_conjDes').set_weights(conjFromWeights);
    n.so('p_stateSim_conjDes').set_weights(conjFromWeights);
    n.so('p_next_conjDes').set_weights(conjToWeights);
end


GOAL = 3;
n.so('n_env').set_vals(e.get_locVec.*4);
n.so('n_goalOut').set_vals(GOAL); n.so('n_goalOut').scaleVals(2);
n.so('n_qOsc').turnOff();

moved=0;
endTime=250; 
for t = dt:dt:endTime
    
    
    %===================
    %=== Exploration ===
    %===================
    if 1
        if t==dt
            subBgFunc=@(vals,in) (any(vals>0)*max((vals-(dt/5)),0));
            n.addNode('n_subBg', numActions)
            n.so('n_subBg').updateFunc=subBgFunc;
            n.so('n_subBg').startLogging();
            n.so('prevMotor1').nodes={n.so('n_subBg')};
        end
        if atTime(10:10:endTime)
            t
            [moved, action]= e.rwalk2();
            n.so('n_env').set_vals(e.get_locVec.*4);
            n.setNode('n_novInhib',1);
            
            
            n.setNode('n_subBg', [1:numActions] == action)
        end
    end
    
    
    %=========================
    %=== Queue Performance ===
    %=========================
    if 0
        if atTime(4) %start queue loading
            n.so('n_qOsc').set_offset(t);
            n.so('n_qOsc').turnOn();
            n.so('n_qOn').set_vals(1);
    %         n.setNode('n_novInhib',1);
        elseif atTime(25) %run queue
            n.so('n_qOn').set_vals(2);
        end

        moved=e.mwalk(n.so('n_bg').vals);
        if moved==1
            n.so('n_env').set_vals(e.get_locVec.*4);
            n.setNode('n_novInhib',1);
        end
    end
%     
    n.update();
  
end

if CLOSE_ALL
    close all;
end

% %State Items
% plotItems = {'state','adjFuture', 'proxPast','next', 'n_novInhib','stateSim', 'prevState1','prevState2', 'goalOut'};
% n.plot(plotItems,50)
% 
% % %Queue Items
% % plotItems={'state','stateSim','n_qOn','next','conjDes','qIn','qStore', 'qOut', 'conjOut', 'n_qOsc'};
% % n.plot(plotItems,50)
% 
%Motor Items
plotItems={'state', 'next', 'conjDes', 'conjObs', 'conjOut', 'motorIn', 'prevMotor1', 'n_bg','prevMotor2', 'n_subBg'};
n.plot(plotItems,50)
% 
plotItems = {'p_conjOut_motorIn'};
n.plot(plotItems, 50);

% 
if 1
% n.plot2({'p_proxPast_next'},50,'proj')

% figure()
% n.projPlot('p_conjDes_qIn','n_qOsc',50)
% 
% figure()
% n.projPlot('p_next_conjDes','n_qOsc',50)
% [x,y]=n.getGateIntervals('p_proxPast_next','n_subOsc',50)
% 
hold off
figure(11)
futureWeights = n.so('p_state_adjFuture').weights;
imagesc(futureWeights);
hold on
title('adjFuture Weights');

hold off
figure(14)
futureWeights = n.so('p_stateSim_adjFuture').weights;
imagesc(futureWeights);
hold on
title('stateSim adjFuture Weights');
% 
% % 
hold off
figure(21)
pastWeights=n.so('p_proxPast_proxPast').weights;
imagesc(pastWeights)
hold on
title('proxPast Weights');


hold off
figure(71)
pastWeights=n.so('p_conjOut_motorIn').weights;
imagesc(pastWeights)
hold on
title('conjOut Motor Weights');


figure(); imagesc(conjActWeights)

% 
% figure(22)
% hold on
% plot(n.so('p_proxPast_proxPast').traceLog);
% n.so('p_proxPast_proxPast').plotGates();
% title('proxPast Traces');
% 
% figure(23)
% hold on
% plot(n.so('proxPast').actLog)
% n.so('p_proxPast_proxPast').plotGates();
% title('proxPast Traces')


% 
% %     
% figure(99)
% imagesc(e.adjacencyMat+eye(e.numStates))
% % 
% figure(12)
% hold on
% plot(n.so('p_state_adjFuture').traceLog)
% title('adjFuture Traces')
end

    
    