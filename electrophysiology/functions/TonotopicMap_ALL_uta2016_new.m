function TonotopicMap_ALL_uta2016(FILEname, fig_flag_mat, flag_parfor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2016/10/28 by Tomoyo
% from Ishizu TonotopicMap_ALL_uta2016(FILEname, MichiganLFP=LFPwaves, SpikeTiming, StmCondMatrix ,fig_flag)
% This program returns TC calcuration result from the data of Michigan electrode.
%   fig_flag = 1;   PSTH
%   fig_flag = 2;   TuneCurve picture
%   fig_flag = 3;   Raser plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tonotopic_map_all_2016(FILEname,flag,analyze_flag)%,available_channel_mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2016/5/26 by Tomoyo
% flag 0[all] 1[Load & Draw Raster] 2[Load & Draw CF] 3[Load & Draw PSTH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

s=['cd TC_',FILEname(16:18)];
eval(s);
name = FILEname(1,9:18);
FreqNum= 18; dBNum= 7; chNum = 96;
Freq = [1.6,2.0,2.5,3.2,4.0,5.0,6.4,8.0,10,13,16,20,25,32,40,50,57,64]; dB = 20:10:80;
%Dur=600;
%Dur_LFP=200;
%CHNUM =96;
%FrqVec = [1.6,2.0,2.5,3.2,4.0,5.0,6.4,8.0,10,13,16,20,25,32,40,50,57,64];
%dBVec  = [80, 70, 60, 50, 40, 30, 20];

%data load
flag=1;
load(FILEname)
if size(find(flag==1))~=0
    %active channel load
    loadNAME=[FILEname,'_activechannel.mat'];
    load(loadNAME)
end
% active_channel(1,96)
% CHdata(1,96).LFP(1,***) .TEMPtimestamps(1,***)
% serial_data(1,10081or10080)
% TRIGwave(1,***about1800000=30min)

%%% make StmCondMatrix
recSCMName = [FILEname(1:18) '_StmCondMatrix.mat'];
if( exist(recSCMName,'file')==0 ),
    disp('make mat file of StmCondMatrix.');
    %            StmCondMatrix = makeStmPTN(serial_data,TRIGwave,preceding_TRIG)% Timestamps = 30kHz
    StmCondMatrix = makeStmPTN_withInput(serial_data,TRIGwave,flag);% Timestamps = 30kHz
    diff(StmCondMatrix(3,61:66))
    disp('It is OK if the third number is larger than others...')
    recSCMName = [FILEname(1:18) '_StmCondMatrix.mat'];
    save( recSCMName, 'StmCondMatrix' );
else
    disp('load StmCondMatrix.');
    load(recSCMName);
end


%%% parameter as to timespan for drawing each figure %%%
PSTHband    = 0:300; % Timespan[ms] from stimulus for drawing PSTH
Raster_band = 0:100; % Raster plot timespan[ms] from stimulus
onset_band  = 5:55;  % Spiking activity timespan[ms], based on stimulus, when neurons are excited in response to the stimulus
spont_band  = 0:4;   % Spontaneous activity timespan[ms], based on stimulus, when neurons are not excited in response to the stimulus
AEPband    = 0:200; % Timespan[ms] from stimulus for drawing AEP
onset_band_long  = 51:100;

%%% make AERwave(7,18,20,100,CHNUM); by Tomoyo %% optional
AERwave=zeros(7,18,20,length(AEPband),chNum);
LOADname = [FILEname(1:18) '_ns2.mat'];
load(LOADname,'analogWaves');
% analogWaves(CHNUM,30min)
if(size(find(fig_flag_mat==0))>0)
    AER_temp=zeros(chNum,length(AEPband),size(StmCondMatrix,2));
    for i=1:size(StmCondMatrix,2)
        AER_temp(:,:,i)=analogWaves(:,StmCondMatrix(3,i)/30:StmCondMatrix(3,i)/30+length(AEPband)-1);
        AER_temp(:,:,i)=AER_temp(:,:,i)-repmat(AER_temp(:,1,i),1,length(AEPband));
    end
    for frq_num=1:18
        for dB_num=1:7
            STMtemp= StmCondMatrix(1,:)==frq_num&StmCondMatrix(2,:)==dB_num;
            AERwave(dB_num,frq_num,:,:,:)=permute(AER_temp(:,:,STMtemp),[3,2,1]);
        end
    end
    clear AER_temp
    recSCMName = [FILEname(1:18) '_AERwave.mat'];
    save( recSCMName, 'AERwave' );
end
%%% spike timing: down sampling from 30kHz to 1kHz %%%
% Probably, sum number of spike are completely preserved by this method.
% Around electrode, firing rate over 1 kHz is not observed.
SpikeTime = zeros(chNum,size(analogWaves,2));
%SpikeTime = zeros(size(MichiganLFP));
for id = 1:size(SpikeTime,1)
    tmp = [];
    try   tmp = round(SpikeTiming{id,1}/30);
    catch
    end
    SpikeTime(id,tmp(tmp>0))=1;
end

%%% trial number in each trial (one freqency/one dB): default => 20 %%%
TrialNum = floor(length(StmCondMatrix)/(FreqNum*dBNum));
disp(['TrialNum =', num2str(TrialNum)]);

%%% arrangement in a sequence of StmCondMatrix according to the Frequency/dB %%%
% result: StmCondMatrix2(1,:)= 1,1,1,1,1....1,1,2,2,2...18 :frequency
%         StmCondMatrix2(2,:)= 1,1,1..2,2,..7,7,1,1......7 :dB
%         StmCondMatrix2(3,:)= ........................... :trigger timing (30kHz sampling)
%
[Disused,idx] = sort(StmCondMatrix(1,:)); tmpStmCondMatrix = StmCondMatrix(:,idx);  tmp  = diff(tmpStmCondMatrix,1,2);
StartId = [1,find(tmp(1,:))+1]; EndId = [find(tmp(1,:)),length(StmCondMatrix)]; dBmatrix = tmpStmCondMatrix(2,:);
idx2 = zeros(1,FreqNum*dBNum*TrialNum);
for n=1:FreqNum
    [dBid, tmpid] = sort(dBmatrix(1,StartId(n):EndId(n)));
    if(length(tmpid) > TrialNum*dBNum);
        tmpstart=[1,find(diff(dBid))+1]; tmpend=tmpstart+TrialNum-1; tmpid2 = zeros(1,TrialNum*dBNum);
        for m = 1:dBNum
            tmpid2(1,TrialNum*(m-1)+1:TrialNum*m) = tmpid(1,tmpstart(1,m):tmpend(1,m));
        end
        idx2(1,TrialNum*dBNum*(n-1)+1:TrialNum*dBNum*n) = tmpid2+TrialNum*dBNum*(n-1);
    else
        idx2(1,TrialNum*dBNum*(n-1)+1:TrialNum*dBNum*n) = tmpid+TrialNum*dBNum*(n-1);
    end
end
StmCondMatrix2 = tmpStmCondMatrix(:,idx2);

% draw figure
for fig_flag=fig_flag_mat
    switch fig_flag,
        case 1, % drawing PSTH
            %%% configuration of each data for drawing figure %%%
            [PSTH] = makePSTHData(SpikeTime, StmCondMatrix2, PSTHband, chNum, TrialNum);
            [subplotnum,plotnum,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID]=PSTHidx(Freq, dB, dBNum, FreqNum);
            SAVEname = [FILEname(1:18) '_PSTH.mat'];
            save(SAVEname,'PSTH');
            
            mkdir PSTH; cd PSTH;
            clc; disp('PSTH drawing.');
            switch flag_parfor
                case 1
                    parfor_progress(chNum); % initializing 'parfor_progress' function
                    parfor ch=1:chNum,
                        recPSTHname_png = [name,'_ch',num2str(ch),'_PSTH.png'];
                        drawPSTH(PSTH,onset_band,ch,dBNum,FreqNum,subplotnum,plotnum,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID)
                        set(gcf,'Name',FILEname);
                        saveas(gcf,recPSTHname_png,'png');
                        parfor_progress; % display the calculation progress
                    end
                    parfor_progress(0); % cleaning up 'parfor_progress' function
                case 0
                    for ch=1:chNum,
                        recPSTHname_png = [name,'_ch',num2str(ch),'_PSTH.png'];
                        drawPSTH(PSTH,onset_band,ch,dBNum,FreqNum,subplotnum,plotnum,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID)
                        set(gcf,'Name',FILEname);
                        saveas(gcf,recPSTHname_png,'png');
                    end
            end
            cd ..
            
        case 2,  % drawing Tune Curve picture
            %%% configuration of each data for drawing figure %%%
            [TCmap, BFmap] = makeTCData(SpikeTime,StmCondMatrix2,onset_band,spont_band,Freq,FreqNum,dBNum,chNum,TrialNum);
            [TCmap_long, BFmap_long] = makeTCData(SpikeTime,StmCondMatrix2,onset_band_long,spont_band,Freq,FreqNum,dBNum,chNum,TrialNum);
            [FRA]=makeFRA(TCmap);
            [FRA_long]=makeFRA(TCmap_long);
            SAVEname = [FILEname(1:18) '_TCmap.mat'];
            save(SAVEname,'TCmap','BFmap','TCmap_long','BFmap_long','FRA','FRA_long');%save('BFmap.mat','BFmap');
            
            mkdir TuneCurve;  cd TuneCurve;
            clc;  disp('TuneCurve figure drawing.');
            switch flag_parfor
                case 1
                    parfor_progress(chNum); % initializing 'parfor_progress' function
                    parfor ch=1:chNum,
                        recTCname_png = [name,'_ch',num2str(ch),'_TC.png'];
%                         TuneCurveDraw(TCmap,ch,Freq,FreqNum)% TuneCurveDraw(TCmap,freqNum,mindB,maxdB,stepdB)
                        TuneCurveDraw_withFRA(TCmap,TCmap_long,FRA,FRA_long,ch,Freq,FreqNum)% TuneCurveDraw(TCmap,freqNum,mindB,maxdB,stepdB)
                        set(gcf,'Name',FILEname);
                        saveas(gcf,recTCname_png,'png');
                        close all
                        parfor_progress; % display the calculation progress
                    end
                    parfor_progress(0); % cleaning up 'parfor_progress' function
                case 0
                    for ch=1:chNum,
                        recTCname_png = [name,'_ch',num2str(ch),'_TC.png'];
%                         TuneCurveDraw(TCmap,ch,Freq,FreqNum)% TuneCurveDraw(TCmap,freqNum,mindB,maxdB,stepdB)
                        TuneCurveDraw_withFRA(TCmap,TCmap_long,FRA,FRA_long,ch,Freq,FreqNum)% TuneCurveDraw(TCmap,freqNum,mindB,maxdB,stepdB)
                        set(gcf,'Name',FILEname);
                        saveas(gcf,recTCname_png,'png');
                        close all
                    end                    
            end
            cd ..
            
        case 3, % drawing Raster plot
            %%% configuration of each data for drawing figure %%%
            [FirstSpikeRaster,Raster,FirstSpikeLatency] = makeRasterData(SpikeTime,StmCondMatrix2,Raster_band,spont_band,FreqNum,dBNum,chNum,TrialNum);
            [plot_id,plotFSR_id,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID]=Raster_idx(FirstSpikeRaster,Raster,Freq,FreqNum,dB,dBNum,TrialNum);
            SAVEname = [FILEname(1:18) '_RasterAndFSL.mat'];
            save(SAVEname,'FirstSpikeRaster','Raster','FirstSpikeLatency');

            mkdir Raster;  cd Raster;
            clc; disp('Raster drawing.');
            switch flag_parfor
                case 1
                    parfor_progress(chNum); % initializing 'parfor_progress' function
                    parfor ch=1:chNum,
                        recRasterplot_png = [name,'_ch',num2str(ch),'_Raster.png'];
                        recRasterplot_fig = [name,'_ch',num2str(ch),'_Raster.fig'];
                        drawRaster(FirstSpikeRaster,Raster,ch,FreqNum,dBNum,TrialNum,plot_id,plotFSR_id,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID);
                        set(gcf,'Name',FILEname);
                        saveas(gcf,recRasterplot_png,'png');
                        savefig(gcf,recRasterplot_fig);
                        close all
                        parfor_progress; % display the calculation progress
                    end
                    parfor_progress(0); % cleaning up 'parfor_progress' function
                case 0
                    for ch=1:chNum,
                        recRasterplot_png = [name,'_ch',num2str(ch),'_Raster.png'];
                        recRasterplot_fig = [name,'_ch',num2str(ch),'_Raster.fig'];
                        drawRaster(FirstSpikeRaster,Raster,ch,FreqNum,dBNum,TrialNum,plot_id,plotFSR_id,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID);
                        set(gcf,'Name',FILEname);
                        saveas(gcf,recRasterplot_png,'png');
                        savefig(gcf,recRasterplot_fig);
                        close all
                    end
            end
            cd ..
    end
end
cd ..
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PSTH] = makePSTHData(SpikeTime, StmCondMatrix2, PSTHband, chNum, TrialNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% index of PSTH timespan spike data %%% PSTHsample*trial*dB*freqency
idxPSTH = repmat(PSTHband,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(PSTHband)]);
idxPSTH = reshape(idxPSTH.',1,[]);

%%% spike data train of each timespan as PSTH, onset timespan %%%
PSTH    = cell(chNum,1); % spike data for PSTH (0-300 ms)
for ch=1:chNum
    tmpPSTHData = reshape(SpikeTime(ch,idxPSTH),length(PSTHband),TrialNum,[]);
    tmpPSTHData = squeeze(sum(tmpPSTHData,2)); % compression in the direction of Trial
    PSTH{ch,1}  = tmpPSTHData; % (0-300ms each trial's hist data) �~ (dBNum*FreqNum)
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TCmap, BFmap] = makeTCData(SpikeTime,StmCondMatrix2,onset_band,spont_band,Freq,FreqNum,dBNum,chNum,TrialNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% index of onset timespan spike data %%% onset_sample*trial*dB*freqency
idx_onset = repmat(onset_band,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(onset_band)]);
idx_onset = reshape(idx_onset.',1,[]);

%%% index of spont timespan spike data %%% spont_band*trial*dB*freqency
idxSpont = repmat(spont_band,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(spont_band)]);
idxSpont = reshape(idxSpont.',1,[]);

%%% firing rate data of onset timespan %%%
onsetSpikeData = zeros(dBNum,FreqNum,chNum);% spike data for onset(5-55  ms)
BFmap = zeros(dBNum,chNum);%
for ch=1:chNum
    tmp_onset = reshape(SpikeTime(ch,idx_onset),length(onset_band),TrialNum, dBNum, FreqNum);% 4-D matrix(datalength,Trial, dB, Freq)
    tmp_onset = squeeze(sum(tmp_onset,2)); % compression in the direction of Trial
    tmp_onset = squeeze(sum(tmp_onset,1)); % compression in the direction of each datalength
    onsetSpikeData(:,:,ch) = tmp_onset / (length(onset_band)*TrialNum); % Averaging to firing rate in the spontaneous duration(5-55ms)
    [M,I]=max(tmp_onset,[],2);
    BFmap(M>0,ch)=Freq(I(M>0));
end

%%% spike data train of each timespan as spontaneous activity %%%
SpontmeanFiring = zeros(dBNum,FreqNum,chNum);% spike data of 0-4 ms after tone stimulus
for ch=1:chNum
    tmp_sponta = reshape(SpikeTime(ch,idxSpont),length(spont_band),TrialNum, dBNum, FreqNum);% 4-D matrix(datalength,Trial, dB, Freq)
    tmp_sponta = squeeze(sum(tmp_sponta,4)); % compression in the direction of FreqNum
    tmp_sponta = squeeze(sum(tmp_sponta,3)); % compression in the direction of dBNum
    tmp_sponta = squeeze(sum(tmp_sponta,2)); % compression in the direction of Trial
    tmp_sponta = squeeze(sum(tmp_sponta,1)); % compression in the direction of each datalength
    spont_firing = tmp_sponta / (dBNum*FreqNum*length(spont_band)*TrialNum); % Averaging to firing rate in the spontaneous duration(0-4ms)
    SpontmeanFiring(:,:,ch) = repmat(spont_firing,dBNum,FreqNum);
end

%%% Tune Curve plot data for CF mapping %%%
TCmap = onsetSpikeData-SpontmeanFiring;% Tune Curve mapping data
TCmap(TCmap<0) = 0 ;


%%% whole datalength firing rate per each channel %%%
%
% MEANband    = 0:599; %[ms]
% meanFiring = zeros(dBNum,FreqNum,chNum);   % spike data of 0-599 ms after tone stimulus
% idxMEAN = repmat(MEANband,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(MEANband)]);
% idxMEAN = reshape(idxMEAN.',1,[]);
% for ch=1:chNum
%     tmp_mean = reshape(SpikeTime(ch,idxMEAN),length(MEANband),TrialNum, dBNum, FreqNum);% 4-D matrix(datalength,Trial, dB, Freq)
%     tmp_mean = squeeze(sum(tmp_mean,2)); % compression in the direction of Trial
%     tmp_mean = squeeze(sum(tmp_mean,1)); % compression in the direction of each datalength
%     meanFiring(:,:,ch) = tmp_mean / (length(MEANband)*TrialNum); % Averaging to firing rate in the spontaneous duration(0-599ms)
% end
%
% TCmap = onsetSpikeData-meanFiring;% Tune Curve mapping data
%
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FirstSpikeRaster,Raster,FirstSpikeLatency] = makeRasterData(SpikeTime,StmCondMatrix2,Raster_band,spont_band,FreqNum,dBNum,chNum,TrialNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% index of timespan 0f 0-40 ms spike histgram for defining first spike %%% FirstSpike_band*trial*dB*freqency
FirstSpike_band = 0:40;
idxFirstSpikeRaster = repmat(FirstSpike_band,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(FirstSpike_band)]);
idxFirstSpikeRaster = reshape(idxFirstSpikeRaster.',1,[]);

%%% index of Raster timespan spike data %%% Raster_band*trial*dB*freqency
idxRaster = repmat(Raster_band,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(Raster_band)]);
idxRaster = reshape(idxRaster.',1,[]);

%%% index of spont timespan spike data %%% spont_band*trial*dB*freqency
idxSpont = repmat(spont_band,[length(StmCondMatrix2),1]) + repmat(round(StmCondMatrix2(3,:)/30).',[1,length(spont_band)]);
idxSpont = reshape(idxSpont.',1,[]);

%%% spike data train of each timespan as Raster %%%
Raster     = cell(dBNum,FreqNum,chNum); % spike data for Rater (0-100 ms)
RasterOnset = cell(chNum,1); % spike data for defining first spike (0-40 ms)
for ch=1:chNum
    tmpRasterData = reshape(SpikeTime(ch,idxRaster),length(Raster_band),[]);
    tmpRasterOnsetData = reshape(SpikeTime(ch,idxFirstSpikeRaster),length(FirstSpike_band),TrialNum,[]);
    tmpRasterOnsetData = squeeze(sum(tmpRasterOnsetData,2)); % compression in the direction of Trial
    tmpRasterCell = mat2cell(tmpRasterData,length(Raster_band),repmat(TrialNum,[1,dBNum*FreqNum]));
    Raster(:,:,ch) = reshape(tmpRasterCell,[dBNum,FreqNum]);
    RasterOnset{ch,1} = tmpRasterOnsetData; % (0-40ms each trial's hist data) �~ (dBNum*FreqNum)
end

%%% spike data train of each timespan as spontaneous activity %%%
sponta = zeros(dBNum,FreqNum,chNum);   % spike data of 0-4 ms after tone stimulus
spikenum  = ones(dBNum,FreqNum,chNum); % spike number neccessary to define first spike latency
for ch=1:chNum
    tmp_sponta = reshape(SpikeTime(ch,idxSpont),length(spont_band),TrialNum, dBNum, FreqNum);% 4-D matrix(datalength,Trial, dB, Freq)
    tmp_sponta = squeeze(sum(tmp_sponta,2)); % compression in the direction of Trial
    tmp_sponta = squeeze(sum(tmp_sponta,1)); % compression in the direction of each datalength
    sponta(:,:,ch) = tmp_sponta / (length(spont_band)*TrialNum); % Averaging to firing rate in the spontaneous duration(0-5ms)
    buff = zeros(dBNum,FreqNum);
    for i=1:TrialNum;
        buff = buff + (sponta(:,:,ch)*20).^(i-1).*exp(sponta(:,:,ch)*(-20))/factorial(i-1);
        pvalue = ones(dBNum,FreqNum) - buff; spikePlus = ones(dBNum,FreqNum);
        spikePlus(pvalue<0.01) = 0; % the first spike depends on the spike number when pvalue goes blow 0.01 at first
        spikenum(:,:,ch) = spikenum(:,:,ch) + spikePlus;
    end
end

%%% calcurating first spike latency (the latency must be within 4-40ms from stimulus)%%%
FirstSpikeLatency  = zeros(dBNum,FreqNum,chNum);
for ch=1:chNum
    tmp = RasterOnset{ch,1}-repmat(reshape(spikenum(:,:,ch),1,[]),[length(FirstSpike_band),1]);
    latency  = zeros(dBNum,FreqNum);
    for i=1:size(tmp,2)
        tmp2 = find(tmp(:,i)>=0);
        try   latency(i) = tmp2(1,1);
        catch
        end
    end
    check1 = zeros(dBNum,FreqNum); check2 = zeros(dBNum,FreqNum);
    check1(latency > 4) = 1; check2(latency < 40) = 1; % the latency only within 4-40 ms is admitted as the CORRECT first spike latency
    FirstSpikeLatency(:,:,ch) = latency.*check1.*check2;
end

%%% making first spike train data (4-40ms from stimulus) %%%
FirstSpikeRaster = cell(dBNum,FreqNum,chNum); % first spike data for Rater (0-40 ms)
for ch=1:chNum
    %%% finding the spike train sfter first spike latency %%%
    tmpcheck = FirstSpikeLatency(:,:,ch);
    check = zeros(length(FirstSpike_band),TrialNum,dBNum*FreqNum);
    for j=1:dBNum*FreqNum
        if(tmpcheck(j)~=0)
            check(tmpcheck(j):end,:,j)=1;
        end
    end
    tmpFirstSpikeRasterData = reshape(SpikeTime(ch,idxFirstSpikeRaster),length(FirstSpike_band),TrialNum,[]);
    tmpFirstSpikeRasterData = tmpFirstSpikeRasterData.*check; % spike count after FirstSpikeLatency
    
    %%% leaving only fist spike %%%
    sp_idx=find(tmpFirstSpikeRasterData);
    [row, col]=find(tmpFirstSpikeRasterData);
    [Uni, idx_col, idUni] = unique(col);
    FirstSpikeTrain = zeros(length(FirstSpike_band),TrialNum,dBNum*FreqNum);
    FirstSpikeTrain(sp_idx(idx_col))=1;
    FirstSpikeTrain = reshape(FirstSpikeTrain,length(FirstSpike_band),[]);
    
    tmpFirstSpikeRasterCell = mat2cell(FirstSpikeTrain,length(FirstSpike_band),repmat(TrialNum,[1,dBNum*FreqNum]));
    FirstSpikeRaster(:,:,ch) = reshape(tmpFirstSpikeRasterCell,[dBNum,FreqNum]);
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawPSTH(PSTH, onset_band, ch, dBNum, FreqNum, subplotnum,plotnum,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHch = PSTH{ch}; % datalength(0-300ms) �~ (dBNum * FreqNum)

%%% onset_band(5-55 ms) spikes are changed red color %%%
Red_start = onset_band(1,1);  %  5 ms
Red_end   = onset_band(1,end);% 55 ms
RedPSTH = [zeros(Red_start,size(PSTHch,2)); PSTHch(Red_start+1:Red_end+1,:);...
    zeros(size(PSTHch,1)-Red_end-1,size(PSTHch,2))];

%%% plotting figure %%%
figure('Position',[200 200 1000 800]); % figure's size and position
for i=1:dBNum*FreqNum,
    subplot(dBNum,FreqNum,subplotnum(i));
    hold on
    line(0:300,PSTHch(:,plotnum(i)).',    'LineWidth',1, 'Color',[0 0 1]);% blue
    line(onset_band, RedPSTH(onset_band+1,plotnum(i)).','LineWidth',2, 'Color',[1 0 0]);% red: onset
    try
        set(gca,'YLim',[0 max(max(PSTHch))],'Fontsize',5) ;
    catch
        set(gca, 'YTick', []);
    end
    set(gca,'XLim',[0 300],'Fontsize',5, 'XTick', [100 300]);
    hold off;
    
    if(dBid(i)),
        text(-150,0.5,num2str(dBstrs(i)),'FontName','Arial','FontSize',12,'Rotation',90);
    end
    if(Freqid(i)),
        text(50,-5,num2str(Freqstrs(i)),'FontName','Arial','FontSize',12);
    end
    if(XtitleID(i)),
        x_str = {'Frequency (kHz)'};
        text(40,-10,x_str,'FontName','Arial','FontSize',16);
    end
    if(YtitleID(i)),
        y_str = {'Intensity (dB SPL)'};
        text(-550,-10,y_str,'FontName','Arial','FontSize',16,'Rotation',90);
    end
end
subplot(dBNum,FreqNum,(floor(FreqNum/2)-1));
titlestring = ['PSTH frq/dB:  ch', num2str(ch)];
text(-100,20,titlestring,'FontName','Arial','FontSize',15);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [subplotnum,plotnum,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID]=PSTHidx(Freq,dB,dBNum,FreqNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot INDEX for drawing the PSTH figure %%%
subplotnum = zeros(dBNum*FreqNum,1); plotnum  = zeros(dBNum*FreqNum,1);
dBid       = zeros(dBNum*FreqNum,1); dBstrs   = zeros(dBNum*FreqNum,1);
Freqid     = zeros(dBNum*FreqNum,1); Freqstrs = zeros(dBNum*FreqNum,1);
XtitleID   = zeros(dBNum*FreqNum,1); YtitleID = zeros(dBNum*FreqNum,1);

t=1;
for i=1:FreqNum,
    for j=1:dBNum,
        subplotnum(t)=(j-1)*FreqNum+i;
        plotnum(t)=(i-1)*dBNum+j;
        
        if i==1,
            dBid(t)=1;
            tdB = flip(dB);
            dBstrs(t) = tdB(j);
        end
        if j == dBNum,
            Freqid(t)=1;
            Freqstrs(t) = Freq(i);
        end
        if j==dBNum && i==floor(FreqNum/2-1),
            XtitleID(t)=1;
        end
        if j==floor(dBNum/2+1) && i==1,
            YtitleID(t)=1;
        end
        t=t+1;
    end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TuneCurveDraw(TCmap,ch,Freq,FreqNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCmapCH= TCmap(:,:,ch);
mindB  = 20;%[dB]
maxdB  = 80;%[dB]
stepdB = 10;%[dB]

hfg = figure('Position',[300 100 1000 800]);
subplot('position',[0.1 0.25 0.4 0.5]);
%%% color plot of TCmap histgram %%%
imagesc(TCmapCH);
sizeofCF = size(TCmapCH); sizeX = sizeofCF(1,2); sizeY = sizeofCF(1,1);
axis image; axis tight; axis([0.5 sizeX+0.5 0.5 sizeY+0.5]);
set(get(gca,'Title'), 'String',['channel No.',int2str(ch)],'FontName', 'caribli', 'FontSize', 13);
set(gca, 'YTickMode','manual');
set(gca, 'XTick', 1:FreqNum); set(gca, 'XTicklabel',Freq);
set(gca, 'YTick', 1:sizeY);   set(gca, 'YTicklabel', maxdB:-stepdB:mindB);
set(get(gca, 'xlabel'), 'String', 'Frequency, kHz');
set(get(gca, 'xlabel'), 'FontSize', 13); set(get(gca, 'xlabel'), 'FontName', 'arial');
set(get(gca, 'ylabel'), 'String', 'Intensity, dBSPL');
set(get(gca, 'ylabel'), 'FontSize', 13); set(get(gca, 'ylabel'), 'FontName', 'arial');
shading flat;

%%% gray scaling %%%
CM = zeros(50,3); % gray scaling by 50 steps
for i = 1 : 50,
    for j = 1 : 3,
        CM(i,j) = abs(1/(50-1) * (50-i));
    end
end
colormap(CM);

%%% TCmap coutour plot %%%
makeContour(TCmapCH,hfg,mindB,maxdB,stepdB,Freq);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TuneCurveDraw_withFRA(TCmap,TCmap_long,FRA,FRA_long,ch,Freq,FreqNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCmapCH= TCmap(:,:,ch);
mapData=zeros(size(TCmap,1),size(TCmap,2),4);
mapData(:,:,1)=TCmap(:,:,ch);
mapData(:,:,2)=FRA(:,:,ch);
mapData(:,:,3)=TCmap_long(:,:,ch);
mapData(:,:,4)=FRA_long(:,:,ch);

mindB  = 20;%[dB]
maxdB  = 80;%[dB]
stepdB = 10;%[dB]

hfg = figure('Position',[300 100 1000 800]);
for i=1:4
    subplot(4,2,(i-1)*2+1);%subplot('position',[0.1 0.25 0.4 0.5]);
    %%% color plot of TCmap histgram %%%
    imagesc(mapData(:,:,i));
    sizeofCF = size(mapData); sizeX = sizeofCF(1,2); sizeY = sizeofCF(1,1);
    axis image; axis tight; axis([0.5 sizeX+0.5 0.5 sizeY+0.5]);
    if(i==1)
        set(get(gca,'Title'), 'String',['channel No.',int2str(ch)],'FontName', 'caribli', 'FontSize', 13);
    end
    set(gca, 'YTickMode','manual');
    set(gca, 'XTick', 1:FreqNum); set(gca, 'XTicklabel',Freq);
    set(gca, 'YTick', 1:sizeY);   set(gca, 'YTicklabel', maxdB:-stepdB:mindB);
    set(get(gca, 'xlabel'), 'String', 'Frequency, kHz');
    set(get(gca, 'xlabel'), 'FontSize', 10); set(get(gca, 'xlabel'), 'FontName', 'arial');
    set(get(gca, 'ylabel'), 'String', 'Intensity, dBSPL');
    set(get(gca, 'ylabel'), 'FontSize', 10); set(get(gca, 'ylabel'), 'FontName', 'arial');
    shading flat;
end

%%% gray scaling %%%
CM = zeros(50,3); % gray scaling by 50 steps
for i = 1 : 50,
    for j = 1 : 3,
        CM(i,j) = abs(1/(50-1) * (50-i));
    end
end
colormap(CM);

%%% TCmap coutour plot %%%
makeContour2(TCmapCH,hfg,mindB,maxdB,stepdB,Freq);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeContour2(TCmapCH,hfg,mindB,maxdB,stepdB,Freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFsize = size(TCmapCH);
TCmapCH = TCmapCH(1:CFsize(1), 1:CFsize(2));
CFsize = size(TCmapCH);
x = CFsize(1,2);%18
y = CFsize(1,1);%5
[X0,Y0] = meshgrid(1:1:x,1:1:y);
[X1,Y1] = meshgrid(1:0.05:x,y:-0.05:1);
CF1 = interp2(X0,Y0,TCmapCH,X1,Y1,'spline');
%%add Tomoyo
maxValue=max(max(CF1));
v=[1:3]*0.25*maxValue;
%% add end
figure(hfg);
subplot(4,2,[2,4,6,8]);%subplot('position',[0.55 0.05 0.4 0.9]);
contourf(CF1,v);% modify Tomoyo
%contourf(CF1,4);
[c,r] = size(CF1);
set(gca, 'XGrid', 'on'); set(gca, 'YGrid', 'on');
set(gca, 'XTick', 1:(r/(x-1)):r); set(gca, 'FontSize', 12); set(gca, 'XTicklabel',Freq);
set(gca, 'YTick', 1:c/(y-1):c);   set(gca, 'FontSize', 12); set(gca, 'YTicklabel',mindB:stepdB:maxdB);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeContour(TCmapCH,hfg,mindB,maxdB,stepdB,Freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFsize = size(TCmapCH);
TCmapCH = TCmapCH(1:CFsize(1), 1:CFsize(2));
CFsize = size(TCmapCH);
x = CFsize(1,2);%18
y = CFsize(1,1);%5
[X0,Y0] = meshgrid(1:1:x,1:1:y);
[X1,Y1] = meshgrid(1:0.05:x,y:-0.05:1);
CF1 = interp2(X0,Y0,TCmapCH,X1,Y1,'spline');
%%add Tomoyo
maxValue=max(max(CF1));
v=[1:3]*0.25*maxValue;
%% add end
figure(hfg);
subplot('position',[0.55 0.05 0.4 0.9]);
contourf(CF1,v);% modify Tomoyo
%contourf(CF1,4);
[c,r] = size(CF1);
set(gca, 'XGrid', 'on'); set(gca, 'YGrid', 'on');
set(gca, 'XTick', 1:(r/(x-1)):r); set(gca, 'FontSize', 12); set(gca, 'XTicklabel',Freq);
set(gca, 'YTick', 1:c/(y-1):c);   set(gca, 'FontSize', 12); set(gca, 'YTicklabel',mindB:stepdB:maxdB);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawRaster(FirstSpikeRaster,Raster,ch,FreqNum,dBNum,TrialNum,plot_id,plotFSR_id,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FirstSpikeRasterCH = cell2mat(FirstSpikeRaster(:,:,ch));
RasterCH           = cell2mat(Raster(:,:,ch));
temp2=1:126;
temp=reshape(fliplr(reshape(temp2,[7,18])),[126,1]);
%save('save','plot_id','plot_id_input','plotFSR_id_input','plotFSR_id');

%%% plotting figure %%%
figure('Position',[350 50 1100 800]);
for i=1:dBNum*FreqNum,
    subplot(FreqNum,dBNum,temp(i))
%     subplot(FreqNum,dBNum,i)
    %%% spike without first spike(blue point) %%%
    [row, col] = find(RasterCH(plot_id{i,1},plot_id{i,2}));
    scatter(row, col,'b','SizeData',2);
    
    %%% first spike(blue point) %%%
    [row_r, col_r] = find(FirstSpikeRasterCH(plotFSR_id{i,1},plotFSR_id{i,2}));
    if(~isempty(row_r)),
        hold on;
        scatter(row_r, col_r,'r','SizeData',2);
        hold off;
    end
    set(gca,'YLim',[0 TrialNum+0.8],'Fontsize',5) ;
    set(gca,'XLim',[0 100],'Fontsize',5);
    
    if(dBid(i)),
        text(40,-20,num2str(dBstrs(i)),'FontName','Arial','FontSize',11);
    end
    if(Freqid(i)),
        text(-50,8,num2str(Freqstrs(i)),'FontName','Arial','FontSize',11);
    end
    if(XtitleID(i)),
        x_str = {'Frequency (kHz)'};
        text(-80,-30,x_str,'FontName','Arial','FontSize',14,'Rotation',90);
    end
%     if(YtitleID(i)),
%         y_str = {'Intensity (dB SPL)'};
%         text(-5,-40,y_str,'FontName','Arial','FontSize',14);
%     end
end
subplot(FreqNum,dBNum,(floor(dBNum/2)+1));
titlestring = ['Raster plot:  ch', num2str(ch)];
text(-40,40,titlestring,'FontName','Arial','FontSize',15);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plot_id,plotFSR_id,dBid,dBstrs,Freqid,Freqstrs,XtitleID,YtitleID]=Raster_idx(FirstSpikeRaster,Raster,Freq,FreqNum,dB,dBNum,TrialNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengthFSR = size( cell2mat(FirstSpikeRaster(:,:,1)) ,1)/dBNum;
lengthR   = size( cell2mat(Raster(:,:,1)) ,1)/dBNum;

%%% plot INDEX for drawing the Raster figure %%%
plot_id  = cell(dBNum*FreqNum,2);  plotFSR_id = cell(dBNum*FreqNum,2);
dBid     = zeros(dBNum*FreqNum,1); dBstrs     = zeros(dBNum*FreqNum,1);
Freqid   = zeros(dBNum*FreqNum,1); Freqstrs   = zeros(dBNum*FreqNum,1);
XtitleID = zeros(dBNum*FreqNum,1); YtitleID   = zeros(dBNum*FreqNum,1);

t=1;
for i=1:FreqNum,
    for j=1:dBNum,
        plot_id{t,1}    = lengthR*(j-1)+1:lengthR*j;
        plot_id{t,2}    = TrialNum*(i-1)+1:TrialNum*i;
        plotFSR_id{t,1} = lengthFSR*(j-1)+1:lengthFSR*j;
        plotFSR_id{t,2} = TrialNum*(i-1)+1:TrialNum*i;
        
        if j == 1,
            Freqid(t)=1;
            Freqstrs(t) = Freq(i);
        end
        if i==FreqNum,
            dBid(t)=1;
            tdB = flip(dB);
            dBstrs(t) = tdB(j);
        end
        if j==1 && i==floor(FreqNum/2+1),
            XtitleID(t)=1;
        end
        if j==floor(dBNum/2+1) && i==FreqNum,
            YtitleID(t)=1;
        end
        t=t+1;
    end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [StmCondMatrix]=makeStmPTN_withInput(serial,STMtimestamps,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2008.12 by otsubo
% 2011.9 kaitei by noda
% 2016.5 kaitei by isoguchi
% STMtimestamps = TRIG waveform
%
% StmCondMatrix = [frqPTN; (1-18)
%                      dB; (1-7)
%         timestamp_ofSTM  (Fs;30kHz)]          3*totalSTM matrix
%                                               timestamp_ofSTM is made from frq inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dBNum  = 7;
ISI=600;

if size(find(flag==2))~=0
    load('stm');
    totalSTM = size(stm,2);
else
    totalSTM = floor(double((length(serial)))/4) ;
end
disp('making STMmatrix')

plot(STMtimestamps(1,1:60*1000));

preceding_TRIG=input('Whether there are any preceding TRIGGER?  (Input) Yes=1, No=0');

%StmCondMatrix = zeros(3,totalSTM+1);%%% 20160407 add last row
StmCondMatrix = zeros(3,totalSTM);%%% 20160407 add last row
if size(find(flag==2))~=0
    StmCondMatrix(1:2,:)=stm;
else
    for k=1:totalSTM,
        StmCondMatrix(1,k) = serial(1,4*k-3);                %Frequency
        TEMPdB = serial(1,4*k-1);
        StmCondMatrix(2,k) = dBNum-TEMPdB+1;
    end
end

%% 2016/5/26 %%
STMtime = zeros(1,totalSTM);
TIMEind = find(STMtimestamps> 1000);
if preceding_TRIG==1
    %%%% for preceding TRIG
    for i=2:length(TIMEind)
        temp1=TIMEind(1,i-1);
        temp2=TIMEind(1,i);
        if temp2-temp1>2000
            STMtime(1,1)=TIMEind(1,i);
            i_now=i;
            break;
        end
    end
else
    %%%% without preceding trigger
    STMtime(1,1) = TIMEind(1,1);
    i_now=2;
end
k = 2;
for i=i_now+1:length(TIMEind),
    temp1 = TIMEind(1,i-1);
    temp2 = TIMEind(1,i);
    if temp2-temp1 > 50,
        STMtime(1,k)=temp2;
        k = k+1;
    end
end

STMtime = STMtime * 30 ;           % 1kS/sec >>> 30kS/sec STMtime=stimulation time

disp(totalSTM);
disp(length(STMtime));

% StmCondMatrix(3,:) = [STMtime(1,[1:totalSTM]),STMtime(1,totalSTM)+ISI*30];
StmCondMatrix(3,:) = STMtime(1,[1:totalSTM]);

disp([StmCondMatrix(:,1:5);StmCondMatrix(3,1:5)/30]);
check_MAT=input('Check the StmCondMatrix!  (Input) OK=1, Wrong=0');
switch check_MAT
    case 0%wrong
        clear StmCondMatrix
        disp('Check the program.')
    case 1%OK
        disp('Start analyze')
end


return;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FRA]=makeFRA(TCmap)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TCmap(dBNum,FreqNum,ChNum)
dBNum=size(TCmap,1);
FreqNum=size(TCmap,2);
ChNum=size(TCmap,3);
BinNum=15;
FRA=zeros(dBNum,FreqNum,ChNum);
Gfilt=fspecial('gaussian',[3,3],0.5);

for ch=1:ChNum
    if(sum(sum(TCmap(:,:,ch)>0)))
        %% gaussian filter 3x3
        Smoothed_TCmap=filter2(Gfilt,TCmap(:,:,ch));
        %% Smoothed_TCmap=imfilter(TCmap,Gfilt,'replicate');
        Smoothed_data=reshape(Smoothed_TCmap,1,[]);
        h=histogram(Smoothed_data,BinNum);
        
        %% calculate Threshold (ver Noda-san)
        x=h.BinEdges;
        y=h.Values;
        m_h=smooth(y,4);
        mm_h=smooth(m_h,3);
        a_mm_h=[0,0,diff(diff(mm_h.'))];
        [max_y,max_y_idx]=max(y);
        smth_idx=max_y_idx;
        [Max_a_mm_h,Max_a_mm_h_idx]=max(a_mm_h(1,smth_idx:end-1));
        Thre_idx=smth_idx+Max_a_mm_h_idx-1;
        Threshold=x(1,Thre_idx);
        
        %% FRA=1/0
        Threshold_data=zeros(1,size(Smoothed_data,2));
        Threshold_data(1,Smoothed_data>Threshold)=1;
        FRA(:,:,ch)=reshape(Threshold_data,dBNum,FreqNum);
    end
end
return;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
