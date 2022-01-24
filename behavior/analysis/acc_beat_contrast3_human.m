%% load music data
cnames = ["s75","s100","s200","s400","click","off"];

load('/Users/yoshiki/research/Tlab/code/sound/trigger/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);

spath = "/Users/yoshiki/research/Tlab/behavior/code/accelerator/software/sound/speed/";

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
cands = ["20210707", "20210708"]; %,"20201208","20201209","20201210","20201211"];
allfiles = [];
for cand = cands
    files = dir(cand+"/*.mat");
    allfile = [];
    for i = 1:length(files)
        if endsWith(files(i).name,".mat")
            allfile = [allfile, files(i)];
        end
    end
    allfiles = [allfiles, allfile];
end
%}

%%

path2 = "/Users/yoshiki/research/Tlab/behavior/human/";
load(path2 + "allfiles.mat")
allfiles([2,6,10]) = [];


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
        monoff = ones(size(ontime,1),4);
        acc2 = zscore(acc, 0, 2);
        acc2 = sqrt(acc2(1,:).^2+acc2(2,:).^2+acc2(3,:).^2);
        jerk2 = diff(acc,1,2)/ dt;
        jerk2 = sqrt(jerk2(1,:).^2+jerk2(2,:).^2+jerk2(3,:).^2);

        for i =1:size(ontime,1)
            [M,Ion] = min(abs(stamp - ontime(i,1)- ontime(i,2)/1000));
            [M,Ioff] = min(abs(stamp - offtime(i,1)- (mlns1(offtime(i,3)+1)-offtime(i,2))/1000));
            monoff(i,2) = Ion;
            %monoff(i,3) = floor(Ion + 60/dt);
            monoff(i,3) = Ioff;
            monoff(i,4) = ontime(i,3)+1;
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
                        tstamps = stamp(st2:en) - stamp(st);
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
lnames2 = cands +"% (" + round((cands/100)*132) + " BPM)";
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


%% plot the max amplitude of acc. beat contrast

lnames2 = cands +"%(" + round((cands/100)*132) + ")";

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.45, 0.6]);
temp = squeeze(values(:,1:5,3));

[p,tbl,stats] = kruskalwallis(temp, [], "off");
disp("Kruskalwallis = " + p)  %2.0e-3
disp("friedman = " + friedman(temp,1,"off"))  %0.0018166
result = multcompare(stats, 'CType','bonferroni');

ps = zeros(length(lnames));
for i = 1:size(result,1)
    ps(result(i,1),result(i,2)) = result(i,end);
end

%{
cand = [1,4;2,4;3,4];
ps = zeros(length(lnames));
for i = 1:length(cand)
    ps(cand(i,1),cand(i,2)) = ranksum(temp(:,cand(i,1)), temp(:,cand(i,2)))*3;
end
%}
disp(ps)

sign = 1*(ps < 0.05) + 1*(ps < 0.01) + 1*(ps < 0.001);
h = boxplot(temp);
hold on
set(h,'LineWidth',1.7)
for i = 1:size(temp,1)
    s1 = scatter((1:size(temp,2))+(rand(1,size(temp,2))-0.5)*0.15, temp(i,:), 130,"ko","filled");
end

ups = [1,1,1];
cnt = 1;
cnt2 = 2.2;

%{
xt = 1:4;
yt = 0.063;
cand = [3,4;2,4;1,4];
for i = 1:length(cand)
    p = ps(cand(i,1),cand(i,2));
    sign = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
    if sign > 0
        cnt = cnt + cnt2*ups(i);
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt*(0.98+cnt/20), '-k', 'LineWidth',1.7)
        for s = 1:sign
            plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/1.5)/6, yt*(0.98+cnt/20) + 0.003, '*k',"markersize",18, 'LineWidth',1.7)
        end
        plot([1;1]*cand(i,1),[0.95+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
        plot([1;1]*cand(i,2),[0.95+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
    else
        cnt = cnt + cnt2*ups(i);
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt*(0.98+cnt/20), '-k', 'LineWidth',1.7)
        text(mean(xt([cand(i,1) cand(i,2)]))-0.1, yt*(0.98+cnt/20)+0.004, "n.s.","Fontsize",25)
        plot([1;1]*cand(i,1),[0.95+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
        plot([1;1]*cand(i,2),[0.95+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
    end
end
%}
lnames_ = lnames;
lnames_(end) = "click";
xticklabels(lnames_)

%xtickangle(20)
if apkup == 1
    ylim([0,0.095])
    yticks((0:10:50)*1e-2)
elseif apkup == 2
    ylim([-0.006,0.7])
    yticks((0:20:70)*1e-2)
end
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;
%xlabel("Play speed(BPM)","Fontsize",32)
%ylabel("Acc. beat contrast  ","Fontsize",32)
%title("Maximum", "Fontsize",40)
saveas(gca, path2 + "summary3/max_acc_beat2"+akinds(apkup)+".fig")
saveas(gca, path2 + "summary3/max_acc_beat2"+akinds(apkup)+".png")

%% circular vector
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.2, 1.2]);
for pkup = 1:5
    subplot(5,1,pkup)
    ex = squeeze(values(:,pkup,1))*2*pi;
    ex2 = squeeze(values(:,pkup,2));
    mn = angle(mean(cos(ex) + 1j*sin(ex),1));
    %ex = mod(ex - mn,2*pi);
    %ex(ex > pi) = 2*pi - ex(ex > pi);
    polarscatter(ex,ex2,60,"xk",'LineWidth',2)
    hold on
    mn2 = mean(ex2.*(cos(ex) + 1j*sin(ex)),1);
    polarscatter(angle(mn2),abs(mn2),60,"or","filled")
    %title(lnames(pkup), "fontsize",17)
    if apkup == 2
        rlim([0,0.4])
        rticks([0,0.2,0.4])
    end
    thetaticks(0:45:360)
    thetaticklabels([0:45:360] + "^o")
    ax = gca;
    ax.RAxis.FontSize = 16;
    ax.ThetaAxis.FontSize = 14;
    ax.RAxisLocation = 290;
    %ax.RAxis.FontWeight = "bold";
    ax.ThetaAxis.FontWeight = "bold";
    ax.GridAlpha = 0.5;
end

saveas(gca, path2 + "summary3/acc_beat_peakpos_a2"+akinds(apkup)+".fig")
saveas(gca, path2 + "summary3/acc_beat_peakpos_a2"+akinds(apkup)+".png")
  
%% inner product test

lnames = ["75%","100%","200%","400%","Click"];
cands = [75,100,200,400, 100];

[ex1,ex2] = pol2cart(squeeze(values(:,1:5,1)).*2*pi, squeeze(values(:,1:5,2)));
ex2 = ex2*1j;
mnp = mean(ex1) + mean(ex2);
ex = zeros(size(ex1));

%
for j = 1:size(ex,1)
    for i = 1:size(ex,1)
        if i == j
            continue
        end
        ex(j,:) = ex(j,:) + (real((ex1(j,:) + ex2(j,:)).*(ex1(i,:) - ex2(i,:))))/(size(ex,1)-1);
    end
end
%}

ex = abs(ex);

%[p,tbl,stats] = anova1(ex, [], "off");
disp("anova1 " + anova1(ex, [], "off"))  %1.1496e-06

[p,tbl,stats] = kruskalwallis(ex, [], "off");
disp("Kruskalwallis = " + p)  %3.5067e-06
disp("friedman = " + friedman(ex,1,"off"))  %4.6694e-06
result = multcompare(stats, 'CType','bonferroni');

%
ps = zeros(length(lnames));
for i = 1:size(result,1)
    ps(result(i,1),result(i,2)) = result(i,end);
end
disp(ps)
%}

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.45, 0.6]);
h = boxplot(ex);
set(h,'LineWidth',1.7)
hold on
for i = 1:size(ex,1)
    scatter((1:5)+(rand(1,5)-0.5)*0.15,ex(i,:),130,"k","filled")
end

cand = [1,2;2,3;3,4;2,4];
ups = [0,0,0,1,1,1,1];
cnt = 1;
cnt2 = 2;
%{
ps = zeros(length(lnames));
for i = 1:length(cand)
    ps(cand(i,1),cand(i,2)) = ranksum(ex(:,cand(i,1)), ex(:,cand(i,2)))*3;
end
disp(ps)
%}

sign = 1*(ps < 0.05) + 1*(ps < 0.01) + 1*(ps < 0.001);

%{
xt = 1:4;
yt = 18e-4;

for i = 1:length(cand)
    p = ps(cand(i,1),cand(i,2));
    sign = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
    if sign > 0
        cnt = cnt + cnt2*ups(i);
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt*(0.98+cnt/20), '-k', 'LineWidth',1.7)
        for s = 1:sign
            plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/1.5)/6,yt*(0.98+cnt/20) + 5e-5, '*k',"markersize",18, 'LineWidth',1.7)
        end
        plot([1;1]*cand(i,1),[0.94+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
        plot([1;1]*cand(i,2),[0.94+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
    end
end
%}

xticklabels(lnames)
yticks((0:100:400)*1e-4)
%ylim([-1e-4,yt*1.25])
ylim([0,0.04])
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;
%xlabel("Play speed(BPM)","Fontsize",32)
%ylabel("$\overrightarrow{r}_i \cdot \overline{\displaystyle\sum_{j \neq i }{\overrightarrow{r_j}}}$","Fontsize",25,'Interpreter','latex')
%ylabel("Beat consensus", "Fontsize",32)
%title("Beat consensus", "Fontsize",40)
saveas(gca, path2 + "summary3/acc_beat_peaktest_a2"+akinds(apkup)+".fig")
saveas(gca, path2 + "summary3/acc_beat_peaktest_a2"+akinds(apkup)+".png")

%%

pkup =4;
I = 31;


temp = zeros(size(allfiles,2),length(lnames),length(gap), length(names));
for j = 1:length(mkinds)
    for k = 1:size(allfiles,2)
        eks = indvdata.("no"+k).(mkinds(j));
        eks_ = eks(:,:,1)./eks(:,:,2);
        temp(k,j,:,:) = (eks_-1)./(eks_+1);
    end
end

anv = ones(length(names),length(gap));
ps = ones(length(names),length(gap),length(mkinds),length(mkinds));

for l = 1:length(gap)
    for g = pkup:pkup
        p = anova1(squeeze(temp(:,:,l,g)),lnames,"off");
        %p = kruskalwallis(squeeze(temp(:,:,l,g)),lnames,"off");
        anv(g,l) = p;
        for j1 = 1:(length(mkinds)-1)
            for j2 = j1:length(mkinds)
                d1 = squeeze(temp(:,j1,l,g));
                d2 = squeeze(temp(:,j2,l,g));
                [p,h] = ranksum(d1,d2);
                %[p,h] = signrank(d1,d2);
                ps(g,l,j1,j2) = p;
            end
        end
    end
end

figure()
imagesc(anv)
disp(min(anv,[],2))

disp([gap(I),anv(pkup,I)])
disp(squeeze(ps(pkup,I,:,:)))

cands = [75,100,200,400,100];
%names = cands +"%(" + round((cands/100)*132) + ")";
name = cands +"%";
name(end) = "Click";


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.45, 0.6]);
h = boxplot(squeeze(temp(:,:,I,pkup)));
set(h,'LineWidth',1.7)
hold on
for i = 1:size(temp,1)
    scatter((1:length(cands))+(rand(1,length(cands))-0.5)*0.15, squeeze(temp(i,:,I,pkup)),130,"k","filled")
end

%{
xt = 1:5;
yt = 0.08;
cand = [1,2;1,3;1,4;2,3;2,4;3,4];
cnt = 2;
cnt2 = 5;
for i = 1:length(cand)
    p = ps(pkup,I,cand(i,1),cand(i,2));
    sign = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
    if sign > 0
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt*(0.98+cnt/20), '-k', 'LineWidth',1.7)
        for s = 1:sign
            plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2)/20, yt*(0.98+cnt/20)+0.003, '*k',"markersize",18, 'LineWidth',1.7)
        end
        plot([1;1]*cand(i,1),[0.9+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
        plot([1;1]*cand(i,2),[0.9+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
        cnt = cnt + cnt2;
    end
end
%}
xticklabels(name)
yticks(-0.3:0.1:0.2)
ylim([-0.3,0.15])
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;

%xlabel("Play speed(BPM)","Fontsize",32)
%ylabel("Acc. beat contrast","Fontsize",32)
title("t = 0 [ms]", "Fontsize",40)
saveas(gca,"./summary3/cont_beat-rnd_a2.fig")
saveas(gca,"./summary3/cont_beat-rnd_a2.png")

%%

%load("null_border.mat")

pkup =4;
I = 31;
gap = -30:1:30;

% the beat contrast at 0 degree and border line (5%)
onbt = temp(:,:,I,pkup);

qualify1 = (onbt - border95(:,2) >= 0) + (onbt - border95(:,1) <= 0);
quantify1 = (onbt - border95(:,2)).*(onbt - border95(:,2) >= 0) + (-onbt + border95(:,1)).*(onbt - border95(:,1) <= 0);
qualify2 = (onbt - border99(:,2) >= 0) + (onbt - border99(:,1) <= 0);
quantify2 = (onbt - border99(:,2)).*(onbt - border99(:,2) >= 0) + (-onbt + border99(:,1)).*(onbt - border99(:,1) <= 0);

disp(sum(qualify1))  %4,7,0,1
disp(sum(qualify2))  %3,5,0,0

% the max beat contrast and border line (5%)
intervals = 43./([0.75,1,2,4]*2.2);
onbt = [];
for i = 1:4
    inter = intervals(i);
    rng = [floor(-inter*3/5), ceil(inter*3/5)] + find(gap==0);
    rng2 = [floor(-inter), ceil(inter)] + find(gap==0);
    data = squeeze(temp(:,i,:,4));
    [M,I] = max(data(:,rng(1):rng(2)),[],2);
    onbt = [onbt,M];
end
%}

qualify3 = (onbt - border95(:,2) >= 0) + (onbt - border95(:,1) <= 0);
quantify3 = (onbt - border95(:,2)).*(onbt - border95(:,2) >= 0) + (-onbt + border95(:,1)).*(onbt - border95(:,1) <= 0);
qualify4 = (onbt - border99(:,2) >= 0) + (onbt - border99(:,1) <= 0);
quantify4 = (onbt - border99(:,2)).*(onbt - border99(:,2) >= 0) + (-onbt + border99(:,1)).*(onbt - border99(:,1) <= 0);

disp(sum(qualify3))  %8,10,7,2
disp(sum(qualify4))  %7,7,2,0




