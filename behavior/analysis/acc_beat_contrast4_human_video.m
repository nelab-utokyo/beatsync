%% load music data
cnames = ["s75","s100","s200","s400","click","off"];

load('/Users/TP/research/Tlab/code/sound/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);

spath = "/Users/TP/research/Tlab/sound/material/";

rpnum = [1,1,2,4,1];
rpnum1 = [1,1,1,4,1];
ratio = [1,0.9946,1,1,1];
mlns = [];

aclick = struct();

for i = 1:5
    if i < 5
        [y,Fs] = audioread(spath + "sonata_" + cnames(i) + ".wav");
        tmp = clspeed.(cnames(i))*1000*ratio(i);
    else
        [y,Fs] = audioread(spath + "click_k448.wav");
        tmp = find(diff(y) > 0.5);
        tmp = tmp(1:2:end)'*1000/96000;
    end
    sound = y(:,1);
    mlns = [mlns, round(length(sound)/Fs*1000)];
    beat = zeros(1, length(allbeats2));
    beat(~isnan(allbeats2)) = tmp(allbeats);
    nans = find(isnan(allbeats2));
    for j = 1:length(nans)-1
        beat(nans(j)) = (beat(nans(j)-1) + beat(nans(j)+1))/2;
    end
    beat(end) = beat(end-1) + median(diff(beat));
    %beat = beat(1:2:end);
    interval = round(length(sound)/Fs*1000) + 50*(i==3);
    if rpnum(i) == 1
        aclick.(cnames(i)) = beat;
    elseif rpnum(i) == 2
        aclick.(cnames(i)) = [beat, beat + interval];
    elseif rpnum(i) == 4
        aclick.(cnames(i)) = [beat, beat + interval, beat + interval*2, beat + interval*3];
    end
end

mlns2 = mlns .* rpnum;

mlns1 = mlns .* rpnum1;

%%
%{
cands = ["20220502","20220513"];
allfiles = [];
for cand = cands
    files = dir(cand+"/App*.mat");
    allfile = [];
    for i = 1:length(files)
        if endsWith(files(i).name,".mat")
            allfile = [allfile, files(i)];
        end
    end
    allfiles = [allfiles, allfile];
end
save("allfiles_video.mat","allfiles")

%}

%%

path2 = "/Users/TP/research/Tlab/behavior/human/";
load(path2 + "allfiles_video.mat")
%allfiles([1,5,11]) = [];

%%
mkinds = ["s75", "s100", "s200", "s400", "click"];

indvdata = struct();
akinds = ["acc","jerk","velocity"];
apkup = 2;  % 1 for acc 2 for jerk 3 for velocity

gap = -30:1:30;
onratio = 0.4;

for k = 1:size(allfiles,2)
    files = allfiles(1,k);
    prep = zeros(length(gap),4,2);
    monamp = struct("s75",prep, "s100",prep, "s200",prep, "s400", prep, "click", prep);
    monamp2 = struct("s75",prep, "s100",prep, "s200",prep, "s400", prep, "click", prep);

    for j = 1:length(files)
        load(files(j).folder + "/" + files(j).name)
        dt = mean(diff(stamp));
        disp(files(j).name)
        disp([k,1/dt])
        dt2 = 1/40;
        monoff = ones(size(ontime,1),5);
        acc2 = zscore(acc, 0, 2);
        acc2 = sqrt(acc2(1,:).^2+acc2(2,:).^2+acc2(3,:).^2);
        jerk2 = diff(acc,1,2)/ dt;
        jerk2 = sqrt(jerk2(1,:).^2+jerk2(2,:).^2+jerk2(3,:).^2);

        for i =1:size(ontime,1)
            [M,Ion] = min(abs(stamp - ontime(i,1)+ ontime(i,2)/1000));
            [M,Ioff] = min(abs(stamp - offtime(i,1)- (mlns1(offtime(i,3)+1)-offtime(i,2))/1000));
            monoff(i,2) = Ion;
            %monoff(i,3) = floor(Ion + 60/dt);
            monoff(i,3) = Ioff;
            monoff(i,4) = ontime(i,3)+1;
            monoff(i,5) = stamp(Ion) - (ontime(i,1) - ontime(i,2)/1000);
            if i < size(ontime,1)
                monoff(i+1, 1) = Ioff + 1;
            end
        end
        for i = 4:4
            for n = 1:size(monoff,1)
                for q = 1:1   %1: correct beat, 2: randomized beat(not used)
                    %
                    onclick = aclick.(mkinds(monoff(n,4)))/1000;
                    oninterval = median(diff(onclick(1:length(beat))))*onratio;
                    st = monoff(n,2);
                    st2 = monoff(n,2)-5;
                    en = monoff(n,3);
                    for g = 1:length(gap)
                        beatime = [onclick-oninterval/2+gap(g)*dt2;onclick+oninterval/2+gap(g)*dt2];
                        tstamps = (stamp(st2:en)+stamp(st2-1:en-1))/2 - stamp(st) + monoff(n,5);
                        tmp = zeros(size(tstamps));
                        for b = 1:length(beatime)
                            tmp((tstamps > beatime(1,b))&(tstamps <= beatime(2,b))) = 1;
                        end
                        if i < 4
                            tmp2 = abs(diff(acc(i,st2-1:en)));
                        else
                            if apkup == 1
                                tmp2 = acc2(st2:en);
                            elseif apkup == 2
                                tmp2 = jerk2(st2-1:en);
                            elseif apkup == 3
                                vel2 = cumsum((acc(:,st2:en) - mean(acc(st2:en),2)).*diff(stamp(st2-1:en)),2)*dt;
                                tmp2 = sqrt(vel2(1,:).^2+vel2(2,:).^2+vel2(3,:).^2);
                            end
                        end     
                        onb = mean(tmp2(tmp==1));
                        offb = mean(tmp2(tmp==0));
                        if q == 1
                            monamp.(mkinds(monoff(n,4)))(g,i,:) = monamp.(mkinds(monoff(n,4)))(g,i,:) + reshape([onb,offb],[1,1,2])/2/length(files); 
                        elseif q > 1
                            monamp2.(mkinds(monoff(n,4)))(g,i,:) = monamp2.(mkinds(monoff(n,4)))(g,i,:) + reshape([onb,offb],[1,1,2])/2/length(files);
                        end
                    end
                end
            end
        end
    end
    indvdata.("no"+k) = monamp;
    indvdata.("no2"+k) = monamp2;
end


%% detecting peaks 

names = ["x", "y", "z", "a"];
cands = [75,100,200,400,100];
lnames = cands + "%";
lnames2 = cands + "% (" + round((cands/100)*132) + " BPM)";
lnames2(end) = "Click (132 BPM)";

temp = zeros(size(allfiles,2),length(lnames),length(gap), length(names));
for j = 1:length(mkinds)
    for k = 1:size(allfiles,2)
        eks = indvdata.("no"+k).(mkinds(j));
        eks_ = eks(:,:,1)./eks(:,:,2);
        temp(k,j,:,:) = (eks_-1)./(eks_+1);
    end
end

if apkup == 1
    lim = [-0.1,0.1];
elseif apkup == 2
    lim = [-0.35,0.35];
elseif apkup == 3
    lim = [-0.05,0.05];
end
clor = ["b","k","r"];
clc = [3,3,3,1,1,3,1,1,3,3];
intervals = 40./([0.75,1,2,4,1]*2.2);

figon = true;
values = zeros(size(temp,1),length(lnames),3);
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 1]);
for pkup = 1:5 %length(lnames)
    subplot(5,1,pkup)
    inter = intervals(pkup);
    datas = [];
    pks = [];
    for i = 1:size(temp,1)
        rng = [floor(-inter*3/5), ceil(inter*3/5)] + find(gap==0);
        rng2 = [floor(-inter), ceil(inter)] + find(gap==0);
        %data = wdenoise(squeeze(temp(i,pkup,:,4)));
        data = squeeze(temp(i,pkup,:,4));
        datas = [datas,data];
        [M,I] = max(data(rng(1):rng(2)));
        [M2,I2] = min(data(rng(1):rng(2)));
        pks = [pks;gap(I+rng(1)-1),M];
        values(i,pkup,1) = mod(I + floor(-inter*3/5) - 1, inter)/inter;
        values(i,pkup,2) = M;
        values(i,pkup,3) = M - M2;
    end
    deg = 45*cands(pkup)/100;
    disp(deg)
    for i = 1:size(temp,1) 
        if values(i,pkup,1)*360 > 90 && values(i,pkup,1)*360 < 270
            cll = clor(1);
        else
            cll = clor(3);
        end
        %cll = clor(3);
        plot(gap,datas(:,i),cll)
        hold on
        scatter(pks(i,1),pks(i,2),60,"kx",'LineWidth',2)
    end
    for i = 0:inter:gap(end)
        plot([1,1]*i, lim, "k")
        plot(-[1,1]*i, lim, "k")
    end
    xticks(gap([7,19,31,43,55]))
    xticklabels(round(gap([7,19,31,43,55])/40,2))
    ylim(lim)
    xlim([gap(1),gap(end)])
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.XAxis.FontSize = 15;
    %if pkup == 4
    %    xlabel("time from the beat [ms]","FontSize",20)
    %end
    %if pkup == 3
    %    ylabel("                        Acc. beat contrast","Fontsize",20)
    %end
    title(lnames2(pkup),"fontsize",21)
end

%saveas(gca, path2 + "summary3/acc_beat_contrast_allr2"+akinds(apkup)+".fig")
%saveas(gca, path2 + "summary3/acc_beat_contrast_allr2"+akinds(apkup)+".png")


%}
%% import video pos data 

mkinds = ["s75", "s100", "s200", "s400", "click"];

datasum = struct();
%available subject: [2,3,4,6,7,8,9,10]
for subno = [2,3,4,6,7,8,9,10]
    
    load(allfiles(subno).folder + "/" + allfiles(subno).name)
    load(allfiles(subno).folder + "/" + replace(allfiles(subno).name, "AppTag","detect_human"))
    
    posdata = struct();
    
    for i = 1:size(ontime,1)
        [M,I1] = min(abs(vstamp - ontime(i,1)));
        [M,I2] = min(abs(vstamp - offtime(i,1)));
        if ~isfield(posdata,mkinds(ontime(i,3)+1))
            %the first round
            posdata.(mkinds(ontime(i,3)+1)) = pos(I1:I2,2:3);
            posdata.(mkinds(ontime(i,3)+1)+"_I") = I1:I2;
            posdata.(mkinds(ontime(i,3)+1)+"_d") = ontime(i,2)/1000;
        else
            %the second round
            posdata.(mkinds(ontime(i,3)+1)+"_2") = pos(I1:I2,2:3);
            posdata.(mkinds(ontime(i,3)+1)+"_2I") = I1:I2;
            posdata.(mkinds(ontime(i,3)+1)+"_2d") = ontime(i,2)/1000;
        end
    end
    
    %
    
    tx = vstamp(posdata.(mkinds(2)+"_I"));
    tx = tx - tx(1) + posdata.(mkinds(2)+"_d");
    [coeff,score,latent] = pca(posdata.s100);
    
    [pkd,lcd,width,prom] = findpeaks(score(:,1),'MinPeakProminence',2);
    [pkd2,lcd2,width2,prom2] = findpeaks(-score(:,1),'MinPeakProminence',2);
    
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.5]);
    
    plot(tx, score(:,1),"linewidth",2)
    hold on
    plot(tx(lcd), pkd,'x','MarkerEdgeColor',[.8 0 .5], "linewidth",1.5)
    plot(tx(lcd2),-pkd2,'x','MarkerEdgeColor',[1 .5 0],"linewidth",1.5)
    for i = 1:132
        plot([1,1]*(aclick.s100(i)/1000),[-max(pkd2)-10,max(pkd)+10],"k")
        hold on
    end
    xlim([0,tx(end)])
    ylim([-max(pkd2)-5,max(pkd)+5])
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.XAxis.FontSize = 15;
    title("100%", "fontsize",25)
    xlabel("time [s]","fontsize",25)
    ylabel("head displacement (PCA)","fontsize",25)
    %saveas(gca, allfiles(subno).folder + "/figure/" + replace(replace(allfiles(subno).name, "AppTag",""),".mat",".png"))
    
    %
    beat = aclick.s100/1000;
    dif = mean(diff(beat));
    
    fig = figure();
    datap = [];
    for d = 1:length(lcd)
        [mn,I] = min(abs(beat - tx(lcd(d))));
        datap = [datap, tx(lcd(d))-beat(I)];
    end
    datap(datap < -dif/2) = [];
    datap(datap > dif/2) = [];
    h1 = histogram(datap, -dif/2:dif/9:dif/2);
    tmp = h1.Values/sum(h1.Values);
    datasum.("up"+subno) = tmp;
    datasum.("up2"+subno) = h1.Values;

    datap2 = [];
    for d = 1:length(lcd2)
        [mn,I2] = min(abs(beat - tx(lcd2(d))));
        datap2 = [datap2, tx(lcd2(d))-beat(I2)];
    end
    datap2(datap2 < -dif/2) = [];
    datap2(datap2 > dif/2) = [];
    
    h2 = histogram(datap2, -dif/2:dif/9:dif/2);
    tmp2 = h2.Values/sum(h2.Values);
    datasum.("low"+subno) = tmp2;
    datasum.("low2"+subno) = h2.Values;
    
    subplot(2,1,1)
    b = bar(1:9,tmp,"hist");
    b.FaceColor = [0.8 0 .5];
    xticks(1:2:9)
    xticklabels(["-π","-π/2","0","π/2","π"])
    ylim([0,ceil((max([tmp,tmp2]))*20)/20])
    ax = gca;
    ax.YAxis.FontSize = 23;
    ax.XAxis.FontSize = 23;
    
    subplot(2,1,2)
    b2 = bar(1:9,tmp2,"hist");
    b2.FaceColor = [1 .5 0];
    ylim([0,ceil((max([tmp,tmp2]))*20)/20])
    xticks(1:2:9)
    xticklabels(["-π","-π/2","0","π/2","π"])
    %histogram(datap2, -dif/2:dif/9:dif/2,"FaceColor",[.8 0 .5])
    
    ax = gca;
    ax.YAxis.FontSize = 23;
    ax.XAxis.FontSize = 23;
    xlabel('Phase');
    %{
    %bar(0.78:1:8.78,tmp,0.44,"FaceColor",[1 .5 0])
    %hold on
    %bar(1.22:1:9.22,tmp2,0.44,"FaceColor",[.8 0 .5])
    %b = bar(tmp2', 1, 'stacked','FaceColor','flat');
    
    xlim([0.45,9.55])
    ylim([0,ceil((max([tmp,tmp2]))*20)/20])
    yticks(0:0.1:0.4)
    %xticks([0,9/4,9/2,9*3/4,9]+0.5)
    xticks([0,8/4,8/2,8*3/4,8]+1)
    xticklabels(["-π","-π/2","0","π/2","π"])
    
    ax = gca;
    ax.YAxis.FontSize = 23;
    ax.XAxis.FontSize = 23;
    xlabel("Phase","Fontsize",26)
    ylabel("Probability","Fontsize",26)
    %}
    
    han=axes(fig,'visible','off'); 
    ax = gca;
    ax.YAxis.FontSize = 25;
    han.YLabel.Visible='on';
    ylabel(han,'Probability');
    
    %saveas(gca, allfiles(subno).folder + "/figure/" + replace(replace(allfiles(subno).name, "AppTag",""),".mat","_hist.png"))
end


%%
fig = figure();
newcolors = {'#E00','#E80','#EE0', '#0B0','#0EE','#00E','#50E','#A0E'};
colororder(newcolors)
tmp = [];
tmp2 = [];
availsub = [2,3,4,6,7,8,9,10];
for subno = availsub
    tmp = [tmp;datasum.("up"+subno)];
    tmp2 = [tmp2;datasum.("low"+subno)];
end
subplot(2,1,1)
b = bar(1:9,tmp/length(availsub),1,"stacked");
%b.FaceColor = [0.8 0 .5];
xticks(1:2:9)
xticklabels(["-π","-π/2","0","π/2","π"])
ylim([0,0.2])
xlim([0,10])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;

subplot(2,1,2)
b2 = bar(1:9,tmp2/length(availsub),1,"stacked");
%b2.FaceColor = [1 .5 0];
ylim([0,0.2])
xlim([0,10])
xticks(1:2:9)
xticklabels(["-π","-π/2","0","π/2","π"])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;
xlabel('Phase');

han=axes(fig,'visible','off'); 
ax = gca;
ax.YAxis.FontSize = 25;
han.YLabel.Visible='on';
ylabel(han,'Probability');

%saveas(gca, "./summary4/video_hist_human.png")

%%
fig = figure();
newcolors = {'#E00','#E80','#EE0', '#0B0','#0EE','#00E','#50E','#A0E'};
colororder(newcolors)
tmp = [];
tmp2 = [];
availsub = [2,3,4,6,7,8,9,10];
interval = linspace(-0.5, 0.5, 10);
interval(end) = [];
for subno = availsub
    %tmp = [tmp;zscore(datasum.("up"+subno)+datasum.("low"+subno))];
    tmp = [tmp;datasum.("up"+subno)+datasum.("low"+subno)];
end

b = bar(interval+0.5/9,tmp,"stacked");

xticks(-0.5:0.25:0.5)
xticklabels(["-π","-π/2","0","π/2","π"])
xlim([-0.5,0.5])
%yticks(-6:3:6)
%ylim([-7,6])
yticks(0:3)
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;

xlabel('Phase', "Fontsize",25);
%ylabel('Zscore',"Fontsize",25);
ylabel('Probability',"Fontsize",25);
%saveas(gca, "./summary4/video_hist_human2.png")

%%


xalpha = reshape(repmat((interval+0.5/9)*pi*2,length(availsub),1),1,[]);
yw = reshape(tmp,1,[])
xalpha = xalpha + pi*(yw < 0);
p = circ_rtest(xalpha, abs(yw))  %0.0080

xalpha = reshape(repmat((interval+0.5/9)*pi*2,length(availsub),1),1,[]);
yw = reshape(tmp/length(availsub),1,[]);
xalpha = xalpha + pi*(yw < 0);
p = circ_rtest(xalpha, abs(yw))  %0.5667
