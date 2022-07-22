%sets up network for small 6-state environment simulation
% Testing modifications to 'find next state' components after April '17 cleanup.
% among other things, changed proxPast learning rule to be oscillatory and
% cleaned up plot method. 

NEW_NETWORK=1;
if NEW_NETWORK
    %clear all
    NEW_NETWORK=1;
else
    n.reset()
end

CLOSE_ALL=1; %closes all figures right before creating new ones (after the run)
SET_WEIGHTS=1; %initialize weights to correct values
ADD_NOISE=0;
rootdir='/data/hammer/space0/model/nz_goals/';
addpath(rootdir);
cd(rootdir);

global dt endTime t 
dt=.05;

START_STATE = 1;

%set up network and environment
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
    
    n.so('proxPast').set_inhibitors({{'n_novInhib',[.8,1,-1],-8},{'n_qOsc',[-.6,-.3,1]},{'n_qOsc',[.6,.3,-1]}});
    
  
    n.addLayer('adjFuture',{numStates,'latInhib',.2,1}); 
    n.so('adjFuture').set_inhibitors({{'n_novInhib',[.7,1,-1]},{'n_qOsc',[-.9,-.3,1]}});
    
    n.addLayer('next',{numStates,'latInhib',1,.1});
    n.so('next').noiseGain=.025;
    n.so('next').bias=-.01;
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
    %n.so('conjObs').set_inhibitors({{'n_novInhib',[.5,.8,0]}});
    n.so('conjObs').set_inhibitors({{'n_novInhib',[.1,1,-1]}});
    
    %=========================
    %== Determining Action ==
    %=========================
    n.addLayer('conjDes', {numStates^2, 'latInhib',1.5,1});
%     n.so('conjDes').set_inhibitors({{'n_novInhib',[.1,1,-1]},{'n_qOsc',[.1, .8, 1]}});
    n.so('conjDes').set_inhibitors({{'n_novInhib',[.1,1,-1]},{'n_qOsc',[.1, 1, 1]},{'n_qOsc',[1, .65, -1]}});

    
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
    n.addOscNode('n_qOsc',1); %1 was working fine; %1.5 was used previously

    %Control signals  
    n.addNode('n_qOn',1); % 0=no cue, 1=storing plan, 2=reading plan
    
    n.addNode('n_noise', 4);
    n.so('n_noise').updateFunc = @(x) min(max(x + dt * (normrnd(0,1,1,4)),0),1);  
    n.addLayer('noise', {numStates, 'dummy', 0, 0});
    n.connectNode('n_noise', 'noise')
    
    %PROJECTIONS
    %(obj,sourceName,targetName,weightConfig,norm,weightFactor,learnType,learnRate,gates)
    %gate format = {'name', [val1, val2]} gate open between val1 & val2
    
    
    %=========================
    %== Choosing the next state ==
    %=========================
    
    n.connect('state','adjFuture','one-to-one',0,2,'adj_future',1,{{'n_qOn',[0 0]}});
    
    n.connectNode('n_goalOut','goalOut'); 
    n.connect('goalOut','proxPast','one-to-one',0,2,'null',0);%, {{'n_subOsc',[-.2,-1]}}); 
   
    n.connect('proxPast','proxPast','one-to-one',2,0,'prox_past',1); 
    n.so('p_proxPast_proxPast').learnGate={{'n_subOsc',[.2,1]}};
    
    n.connect('proxPast', 'next', 'one-to-one', 0,2,'null',0, {{'n_subOsc',[-1,-.2,1]}})%,{'n_qOn',[0 0]}});
    
    
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
    n.connect('conjObs' ,'conjOut', 'impose2', 0, 1, 'null', 0, {{'n_subOsc', [.2, 1]},{'n_qOn', [0,0]}})
    %=========================
    %== Determining Action ==
    %=========================
    
     n.connect('state','conjDes','unif',0,1,'null',0,{{'n_qOn',[0,0]}});    
   

    n.connect('next','conjDes','unif',0,1,'null',0); 
    
    
    n.connect('conjDes','conjOut','impose2',0,1,'null',0, {{'n_subOsc', [-.2,-1]}, {'n_qOn', [0,0]}});
    n.connect('conjDes','conjDes','impose',0,2,'null',0);
    n.connect('conjOut','conjOut','impose',0,2,'null',0);
    n.connect('conjOut', 'motorIn','null',0,0,'oscHebb',0, {{'n_subOsc',[-.2,-1]}}); %learning rate was 20
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
    
    n.connect('qOut','conjOut','one-to-one',0,3,'null',1)%,{{'n_qOsc', [.3,-.3, -1]}});
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
    
    
    n.so('p_conjOut_motorIn').set_weights(conjActWeights);
%     n.so('p_conjOut_motorIn').learnRate=0;
    
    %remove weights onto units representing a "transition" from a state to
    %itself
    conjFromWeights(:,[1, 6, 11, 16]) = 0;
    conjToWeights(:,[1, 6, 11, 16]) = 0;

    n.so('p_state_conjDes').set_weights(conjFromWeights);
    n.so('p_stateSim_conjDes').set_weights(conjFromWeights);
    n.so('p_next_conjDes').set_weights(conjToWeights);
    
end


%random debugging
n.disconnect('conjDes', 'conjOut');
% n.disconnect('conjObs', 'conjOut');


if ADD_NOISE
    noiseVal = .001;
    n.so('state').noiseGain=noiseVal;
    n.so('goalOut').noiseGain=noiseVal;
    n.so('proxPast').noiseGain=.001;
    n.so('adjFuture').noiseGain=noiseVal;
    n.so('prevState1').noiseGain=noiseVal;
    %n.so('prevState2').noiseGain=noiseVal;
    n.so('prevMotor1').noiseGain=noiseVal;
    %n.so('prevMotor2').noiseGain=noiseVal;
    n.so('conjObs').noiseGain=noiseVal;
    n.so('conjDes').noiseGain=noiseVal;
    %n.so('conjOut').noiseGain=noiseVal;
    n.so('motorIn').noiseGain=noiseVal;
    n.so('motorOut').noiseGain=noiseVal;
    n.so('stateSim').noiseGain=noiseVal;
    n.so('qIn').noiseGain=noiseVal;
    n.so('qStore').noiseGain=noiseVal;
    %n.so('qOut').noiseGain=noiseVal;
end



    
    