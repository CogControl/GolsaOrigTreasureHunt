function value = mgRunTrial(net,environment,startState,nextState,finalState) 
    %net = established network
    %environment = state space
    %start = starting state
    %goal1 = next goal
    %goal2 = final goal
    %delay1 = delay after text screen, before response (delay between qLoad
    %and qRead)
    %delay2 = delay between action and outcome (update n_env) 
    %delay3 = ITI (not important?)
    
    %how to ensure goal1 is loaded properly when ultimately heading for
    %goal2?
    global t dt endTime
    n = net;
    e = environment;
    
    n.reset()
    
    %initialize current state and immediate goal
    e.set_loc(startState);
    n.so('n_env').set_vals(e.get_locVec.*4);
    n.so('n_goalOut').set_vals(nextState); n.so('n_goalOut').scaleVals(2);
    n.so('n_qOsc').turnOff();
    n.so('n_qOn').set_vals(0);
    n.so('n_noise').set_vals(ones(1,4) * .5);

    
    close all
    
    moved=0;
    endTime=22;
    endTime = 30;
    for t = dt:dt:endTime

        %=========================
        %=== Queue Performance ===
        %=========================
        if 1
            if atTime(2-dt) %start queue loading
                n.so('n_qOsc').set_offset(t);
                n.so('n_qOsc').turnOn();

        %         n.setNode('n_novInhib',1);
            elseif atTime(4)
                n.so('n_qOn').set_vals(1);

            elseif atTime(7)  %10
               n.so('n_goalOut').set_vals(finalState); n.so('n_goalOut').scaleVals(2);

            elseif atTime(15) %20 %18
               n.so('n_qOn').set_vals(2);
               
           elseif atTime(25)
               n.so('n_qOn').set_vals(2);
            end


            %n.so('n_goalOut').set_vals(GOAL); n.so('n_goalOut').scaleVals(2);
            moved=e.mwalk(n.so('n_bg').vals);
            if moved==1
                n.so('n_env').set_vals(e.get_locVec.*4);
                n.setNode('n_novInhib',1);
            end

        end   
        n.update();
  
    end
    
    
    
%     if e.get_loc() ~= nextState

%         plotTime = 10;
%         %State Items
%         plotItems = {'state','adjFuture', 'proxPast','next', 'n_novInhib','stateSim', 'prevState1','prevState2', 'goalOut'};
%         n.plot(plotItems, plotTime)
%         % 
%         %Queue Items
%         plotItems={'state','stateSim','n_qOn','next','qIn','conjDes','qStore', 'qOut', 'conjOut', 'n_qOsc'};
%         n.plot(plotItems, plotTime)
% 
%         %Motor Items
%         plotItems={'state', 'next', 'conjDes', 'conjObs', 'conjOut', 'motorIn', 'prevMotor1', 'n_bg','prevMotor2', 'n_subOsc'};
%         n.plot(plotItems, plotTime)
        
       
%         error()
%     end
    
% %     Presentation Items
%     figure(44)
%     plotItems = {'state','adjFuture', 'proxPast','next', 'stateSim', 'goalOut', 'qIn', 'qStore', 'qOut', 'n_noise', 'motorIn', 'motorOut', 'noise'};
%     n.plot(plotItems,50)
%     
    value = e.get_loc();

end