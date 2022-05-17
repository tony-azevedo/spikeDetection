% [data,si,h]= abfload('E:\Data\SydneyData\20106001.abf');

a = load('C:\Users\tony\Downloads\SampleForSpikeSorting\SampleForSpikeSorting\021622_Exp4A_bout1_forSpikeSorting.mat')

data = -a.VoltageTrace;

vars_initial.fs = 1000;
vars_initial.spikeTemplateWidth = round(0.005*vars_initial.fs)+1;
vars_initial.len = length(data);
vars_initial.hp_cutoff = 400;
vars_initial.lp_cutoff = 209.2370;
vars_initial.diff = 1;
vars_initial.peak_threshold = 3.7670e-04;
vars_initial.Distance_threshold= 11.0418;
vars_initial.spikeTemplate = [];
vars_initial.likelyiflpntpeak = 33;
vars_initial.polarity = -1;
vars_initial.Amplitude_threshold = 0.7104;

vars_initial.lastfilename = '021622_Exp4A_bout1_forSpikeSorting'

% from a first attempt at spike sorting, gets you started. Comment it out
% to start afresh.
vars_initial.spikeTemplate = ...
     [4.36876936769883e-07,-9.94181993119135e-07,-2.43445885688151e-06,-4.03174272269512e-06,-5.58208662190906e-06,-6.54234827967968e-06,-6.35711529960029e-06,-4.77835347960607e-06,-2.31390590895983e-06,8.67030391281994e-07,4.31874251500089e-06,7.39407239588119e-06,1.00232839956894e-05,1.16176633186079e-05,1.20422905300695e-05,1.18838217366449e-05,1.15545539427135e-05,1.12820488705590e-05,1.05515065544164e-05,8.80410110785578e-06,5.81255852621189e-06,1.36834965101920e-06,-4.90361332999254e-06,-1.42815899705199e-05,-2.85277358892437e-05,-4.88227825530086e-05,-7.52295295409007e-05,-0.000107264931030865,-0.000143995482844807,-0.000183341562111822,-0.000222018994445691,-0.000255988867008909,-0.000281662195580804,-0.000296160384455022,-0.000297148207338936,-0.000283658321489409,-0.000255807131296038,-0.000213932675457374,-0.000159372861486759,-9.41830506323997e-05,-2.04577868664350e-05,5.84923664987400e-05,0.000139352492699211,0.000219550326992173,0.000295572232710245,0.000364563610556225,0.000424814817470689,0.000473729983399184,0.000509700903928211,0.000532583842092633,0.000541888725896611,0.000538107568444536,0.000522631113568304,0.000496811042750134,0.000462407387349446,0.000420935975576188,0.000374124896568031,0.000324247675268120,0.000272903628520082,0.000221471969755249,0.000171412805713445,0.000123366345268516,7.80001536232868e-05,3.61558587565190e-05,-2.28361752500970e-06,-3.75833934658188e-05,-6.96482741228668e-05,-9.85794178182276e-05,-0.000124328516227921,-0.000146631793478390,-0.000165599221037463,-0.000181569213291098,-0.000194650815823694,-0.000204839948495081,-0.000212407272449602,-0.000218045251806509,-0.000222397919786858,-0.000225624670732958,-0.000227772361978792,-0.000229035939495294,-0.000229726245757243,-0.000230073439078513,-0.000229962164789776,-0.000229030028193859,-0.000226768364505393,-0.000222786942737766,-0.000216925175085088,-0.000209301272237913,-0.000200237394382765,-0.000190012838950483,-0.000179106797994864,-0.000168122612035723,-0.000157636183291767,-0.000147876437778514,-0.000138497289898876,-0.000129196361322122,-0.000119739454792159,-0.000109894399626875,-9.96370160100537e-05,-8.94538804692800e-05,-8.00733559204383e-05];


% split data into 20s chunks
% twntysec = 20*vars_initial.fs;
% ttlchnks = floor(length(data)/twntysec);
% data = [data ; zeros(twntysec-(length(data)-(ttlchnks*twntysec)) ,1)];
% chcnkdata = reshape(data,twntysec,[]);

%% Play around, find a nice chunk
% c = 7
% figure
% plot(chcnkdata(:,c))

%% With a chunk you like, find the paramaters for search
[trial,vars_skeleton] = spikeDetection_forABF(data,vars_initial);

%vars_nice = getacqpref('FlyAnalysis',['Spike_params_voltage_fs' num2str(vars_initial.fs)]);


%%
for c = 1:size(chcnkdata,2)

    vars_nice.lastfilename = ['E:\Data\SydneyData\20106001_chnck_' num2str(c) '.mat'];
    [trial,vars_skeleton] = spikeDetection_forABF(chcnkdata(:,c),vars_nice);
    vars_nice = getacqpref('FlyAnalysis',['Spike_params_voltage_fs' num2str(vars_initial.fs)]);

end

% vars_initial.spikeTemplate = vars_skeleton.spikeTemplate;


