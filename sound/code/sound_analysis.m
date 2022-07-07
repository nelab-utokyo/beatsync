%% load music data
cnames = ["s75","s100","s200","s300","s400"];

load('/Users/TP/research/Tlab/code/sound/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);
nonb = 1:489;
nonb(allbeats) = [];

spath = "/Users/TP/research/Tlab/sound/material/";

rpnum = [1,1,2,1,4];
rpnum1 = [1,1,1,1,4];
ratio = [0.9997,0.9946,1,1,1];
%ratio = [1,1,1,1];
mlns = [];

aclick = struct();

for i = 1:5
    [y,Fs] = audioread(spath + "sonata_" + cnames(i) + ".wav");
    sound = y(:,1);
    mlns = [mlns, round(length(sound)/Fs*1000)];
    tmp = clspeed.(cnames(i))*1000*ratio(i);
    beat = zeros(1, length(allbeats2));
    beat(~isnan(allbeats2)) = tmp(allbeats);
    nbeat = tmp(nonb);
    nans = find(isnan(allbeats2));
    for j = 1:length(nans)-1
        beat(nans(j)) = (beat(nans(j)-1) + beat(nans(j)+1))/2;
    end
    beat(end) = beat(end-1) + median(diff(beat));
    interval = round(length(sound)/Fs*1000) + 50*(i==3);
    if rpnum(i) == 1
        aclick.(cnames(i)) = beat;
    elseif rpnum(i) == 2
        aclick.(cnames(i)) = [beat, beat + interval];
    elseif rpnum(i) == 4
        aclick.(cnames(i)) = [beat, beat + interval, beat + interval*2, beat + interval*3];
    end
    aclick.("n"+cnames(i)) = nbeat;
end

mlns2 = mlns .* rpnum;

mlns1 = mlns .* rpnum1;


%%
meanyda = [];
alltmgs = struct();


for snum = 1:5
    [yda,Fs] = audioread(spath + "sonata_" + cnames(snum) + ".wav");
    %x = (1:length(y))/Fs;
    xda = (1:length(yda));
    
    spd = [0.75,1,2,3,4];
    if snum == 1
        ymean = mean(yda(:,1).*yda(:,1));
    end    
    yda(:,1) = yda(:,1)*sqrt(ymean./mean(yda(:,1).*yda(:,1)));
    yda2 = yda(:,1).*yda(:,1);
    [yupper,ylower] = envelope(yda(:,1).*yda(:,1),3000/spd(snum),"peak");

    meanyda = [meanyda, mean(yupper(1:end-Fs,1))];
    ddif = [0,0,0,0,0];
    %timing_ = clspeed.(cnames(2));
    timing_ = round(aclick.(cnames(snum))(1:132)/1000 * Fs)+ddif(snum);
    ntiming_ = round(aclick.("n"+cnames(snum))/1000 * Fs)+ddif(snum);
    timing = [];
    excld = [];
    for i = 1:length(timing_)
        [Mx,Ix] = max(yupper(timing_(i):timing_(i)+Fs/10/spd(snum)));
        %Mx = mean(yda2(timing_(i):timing_(i)+Fs/10));
        %Mx = mean(yda2(timing_(i):timing_(i)+Fs/10));
        Ix = Ix + timing_(i) - 1;
        timing = [timing, Mx];
        %timing = [timing, Mx - yupper(max([1,Ix-0.01*Fs]))];
        excld = [excld, Ix];
    end
    nonbmean = [];
    for i = 1:length(ntiming_)
        [Mx,Ix] = max(yupper(ntiming_(i):ntiming_(i)+Fs/10/spd(snum)));
        %Mx = mean(yda2(ntiming_(i):ntiming_(i)+Fs/10));
        Ix = Ix + ntiming_(i) - 1;
        if sum(excld == Ix) >= 0
            nonbmean = [nonbmean, Mx];
            %nonbmean = [nonbmean, Mx - yupper(max([1,Ix-0.01*Fs]))];
        end
    end
    alltmgs.(cnames(snum) + "x") = timing_/Fs;
    alltmgs.(cnames(snum) + "y") = timing;
    alltmgs.(cnames(snum) + "n") = nonbmean;

    %{
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.25]);
    plot(xda/Fs,yda(:,1),"color",[0.3,0.3,0.3])
   
    xlim([0,60])
    ylim([-1,1])
    xticklabels([])
    yticklabels([])
    axis off
    %[h,icons] = legend("raw waveform","Fontsize",28);
    %icons2 = findobj(icons,'Type','line');
    %set(icons2,'linewidth',3);
    
    %saveas(gca,"./summary/raw_wav.png")

    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.5]);

    plot(xda/Fs,yda(:,1).*yda(:,1),"k")
    hold on
    
    plot(xda/Fs,yupper,"g")
    %scatter(timing_(allbeats),sqrt(mastrg(allbeats,2))/20,"xk")
    plot(timing_/Fs,timing,".r","markersize",15,"linewidth",1.5)
    
    %xlim([0,60])
    ylim([0,0.4])
    yticks(0:0.2:1)
    ax = gca;
    ax.YAxis.FontSize = 28;
    ax.XAxis.FontSize = 28;
    %ylabel("Squared waveform","Fontsize",23)
    %xlabel("Time [s]","Fontsize",28)
    %[h,icons] = legend(["squared waveform","beat"],"Fontsize",28,"Location","north");
    %icons1 = findobj(icons,'Type','line');
    
    %set(icons1,"markersize",15);
    %set(icons1,'linewidth',3);
    box off
    
    %saveas(gca,"./summary/detect_peak.png")
    %}
end

%%

sound4pos = zeros(5,5); %5 speeds 4 positions

for snum = 1:5
    tmp = [mean(alltmgs.(cnames(snum)+"y")(1:4:end)),mean(alltmgs.(cnames(snum)+"y")(2:4:end)),...
        mean(alltmgs.(cnames(snum)+"y")(3:4:end)),mean(alltmgs.(cnames(snum)+"y")(4:4:end))];
    sound4pos(snum,1:4) = tmp;
    sound4pos(snum,5) = mean(alltmgs.(cnames(snum)+"n"));
end

cands = [75,100,200,300,400];
names = cands +"%(" + round((cands/100)*132) + ")";

clors = [0.8500,0.3250,0.0980;0,0.4470,0.7410;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560;0.5,0.5,0.5];

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);

b = bar(sound4pos);

for p = 1:5
    b(p).FaceColor = 'flat';
    b(p).CData = [clors(p,:);clors(p,:);clors(p,:);clors(p,:);clors(p,:)];
end

xticks(1:length(cands))
xticklabels(cands + "%")
%yticks(-0.4:0.1:1)
%ylim([70,75])
[h,icons] = legend(["pos.1","pos.2","pos.3","pos.4","non beat"], "location","north east","fontsize", 23);
xlim([0.5,length(names)+0.5])
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;

ylabel("mean amplitude")
saveas(gca,"./summary/dB_4pos.png")


