cd('/data/hammer/space0/model/nz_goals/tests/mgTask/')

clear all

mgNet;
mgRunTrial(n, e, 1, 2, 2);


plotItems = {'state', 'adjFuture', 'proxPast', 'next', 'stateSim', 'goalOut',...
    'qIn', 'qStore', 'qOut', 'n_bg', 'motorIn', 'motorOut', 'conjDes',...
    'conjObs', 'conjOut'};
times = [.05,30];
fileName = '';

%root: '/data/hammer/space0/model/nz_goals/'
saveDir = 'dissertation/mgTask/trial2/';

n.exportPlotData(plotItems,times,fileName, saveDir);