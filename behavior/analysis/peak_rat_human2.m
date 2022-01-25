%% load music data
cnames = ["s75","s100","s200","s400","click","off"];

load('../../electrophysiology/sound/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);

spath = "../software/sound/speed/";

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


%% calculate peak acc of each beat (rat)

apkup = 1;  %choose the calculation method below
akinds = ["acc","Z(acc)","jerk","Z(jerk)"];

mkinds = ["s75", "s100", "s200", "s400"];

path1 = "/Users/yoshiki/research/Tlab/behavior/acceleration/";
load(path1 + "dno.mat")
allfiles = [dno1,dno2,dno3,dno4,dno5,dno6,dno7,dno8,dno9,dno11];


beatfeatr = zeros(size(allfiles,2),length(mkinds),132);  %collecting the acc value on each beat
maxbeatr = zeros(size(allfiles,2),length(mkinds),132,2); %collecting the maximum acc between beats

certify = zeros(size(allfiles,2),length(mkinds));

for k = 1:size(allfiles,2)
    files = allfiles(1,k);
    for j = 1:length(files)
        load(files(j).folder + "/" + files(j).name)
        dt0 = mean(diff(stamp));
        disp(files(j).name)
        disp([dt0,var(diff(stamp))]*1000)
        monoff = ones(size(ontime,1),4);
        for i =1:size(ontime,1)
            [M,Ion] = min(abs(stamp - ontime(i,1)- ontime(i,2)/1000));
            [M,Ioff] = min(abs(stamp - offtime(i,1)- (mlns1(offtime(i,3)+1)-offtime(i,2))/1000));
            monoff(i,2) = Ion;
            monoff(i,3) = Ioff;
            monoff(i,4) = ontime(i,3)+1;
            if i < size(ontime,1)
                monoff(i+1, 1) = Ioff + 1;
            end
        end
        if apkup == 1
            tmp2_ = acc;
        elseif apkup == 2
            tmp2_ = zscore(acc,0,2);
        elseif apkup == 3
            tmp2_ = diff(acc,1,2)/dt0;
        elseif apkup == 4
            tmp2_ = zscore(diff(acc,1,2)/dt0,0,2);
        end
        tmp2 = sqrt(tmp2_(1,:).^2+tmp2_(2,:).^2+tmp2_(3,:).^2);
        for i = 1:size(ontime,1)
            st1 = monoff(i,2);
            Xst = stamp(st1);
            mpkup = monoff(i,4);
            if mpkup > 4
                continue
            end
            certify(k,mpkup) = certify(k,mpkup) + var(tmp2(st1: monoff(i,3)));
            dt = median(diff(aclick.(mkinds(mpkup))))/1000;
            for i2 = 1:132
                if rem(i2,132) == 0
                    rng = [Xst+aclick.(mkinds(mpkup))(i2)/1000,Xst+aclick.(mkinds(mpkup))(i2)/1000 + dt];
                else
                    rng = [Xst+aclick.(mkinds(mpkup))(i2)/1000,Xst+aclick.(mkinds(mpkup))(i2+1)/1000];
                end
                tmp = (stamp >= rng(1)) & (stamp < rng(2));
                [~,idt] = min(abs(stamp - rng(1)));
                beatfeatr(k,mpkup,rem(i2,132)+1) = beatfeatr(k,mpkup,rem(i2,132)+1) + tmp2(idt)/2;
                if ~isempty(tmp2(tmp))
                    maxbeatr(k,mpkup,rem(i2,132)+1,1) = maxbeatr(k,mpkup,rem(i2,132)+1,1) + max(tmp2(tmp));
                    maxbeatr(k,mpkup,rem(i2,132)+1,2) = maxbeatr(k,mpkup,rem(i2,132)+1,2) + 1;
                end
            end
        end 
    end
end

maxbeatr(:,:,:,1) = maxbeatr(:,:,:,1)./maxbeatr(:,:,:,2);

thresh = nanmean(maxbeatr(:,:,:,1),3) + 2.5*nanstd(maxbeatr(:,:,:,1),0,3);
maxbeatr(maxbeatr(:,:,:,1) - thresh > 0) = NaN;
ratbeat_ = squeeze(nanmean(maxbeatr(:,:,:,1)));
ratbeat = (ratbeat_ - nanmean(ratbeat_,2))./nanstd(ratbeat_,0,2);


%% calculate peak acc of each beat (human)

path2 = "/Users/yoshiki/research/Tlab/behavior/human/";
load(path2 + "allfiles.mat")
allfiles([2,6,10]) = [];

mkinds = ["s75", "s100", "s200", "s400"];
akinds = ["acc","Z(acc)","jerk","Z(jerk)"];

beatfeat = zeros(size(allfiles,2),length(mkinds),132);
maxbeat = zeros(size(allfiles,2),length(mkinds),132,2);

for k = 1:size(allfiles,2)
    files = allfiles(1,k);
    for j = 1:length(files)
        load(files(j).folder + "/" + files(j).name)
        dt0 = mean(diff(stamp));
        disp(files(j).name)
        disp([dt0,var(diff(stamp))]*1000)
        monoff = ones(size(ontime,1),4);
        for i =1:size(ontime,1)
            [M,Ion] = min(abs(stamp - ontime(i,1)- ontime(i,2)/1000));
            [M,Ioff] = min(abs(stamp - offtime(i,1)- (mlns1(offtime(i,3)+1)-offtime(i,2))/1000));
            monoff(i,2) = Ion;
            monoff(i,3) = Ioff;
            monoff(i,4) = ontime(i,3)+1;
            if i < size(ontime,1)
                monoff(i+1, 1) = Ioff + 1;
            end
        end
        if apkup == 1
            tmp2_ = acc;
        elseif apkup == 2
            tmp2_ = zscore(acc,0,2);
        elseif apkup == 3
            tmp2_ = diff(acc,1,2)/dt0;
        elseif apkup == 4
            tmp2_ = zscore(diff(acc,1,2)/dt0,0,2);
        end
        tmp2 = sqrt(tmp2_(1,:).^2+tmp2_(2,:).^2+tmp2_(3,:).^2);
        for i = 1:size(ontime,1)
            st1 = monoff(i,2);
            Xst = stamp(st1);
            mpkup = monoff(i,4);
            if mpkup > 4
                continue
            end
            dt = median(diff(aclick.(mkinds(mpkup))))/1000;
            for i2 = 1:132
                if rem(i2,132) == 0
                    rng = [Xst+aclick.(mkinds(mpkup))(i2)/1000,Xst+aclick.(mkinds(mpkup))(i2)/1000 + dt];
                else
                    rng = [Xst+aclick.(mkinds(mpkup))(i2)/1000,Xst+aclick.(mkinds(mpkup))(i2+1)/1000];
                end
                tmp = (stamp >= rng(1)) & (stamp < rng(2));
                [~,idt] = min(abs(stamp - rng(1)));
                beatfeat(k,mpkup,rem(i2,132)+1) = beatfeat(k,mpkup,rem(i2,132)+1) + tmp2(idt)/2;
                if ~isempty(tmp2(tmp))
                    maxbeat(k,mpkup,rem(i2,132)+1,1) = maxbeat(k,mpkup,rem(i2,132)+1,1) + max(tmp2(tmp));
                    maxbeat(k,mpkup,rem(i2,132)+1,2) = maxbeat(k,mpkup,rem(i2,132)+1,2) + 1;
                end
            end
        end 
    end
end

maxbeat(:,:,:,1) = maxbeat(:,:,:,1)./maxbeat(:,:,:,2);
humanbeat_ = squeeze(nanmean(maxbeat(:,:,:,1)));
humanbeat = zscore(humanbeat_,0,2);


%% plot compared figure of human and rat (each beat)

lnames = ["75%","100%","200%","400%"];

[yda,Fs] = audioread(spath + "sonata_" + cnames(2) + ".wav");
xda = (1:length(yda));

for i = 1:length(lnames)
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.4]);
    if i == 1
        axes('Units','normalized','Position',[0.105 0.75 0.805 0.13]);
        plot(xda/Fs,yda(:,1)/0.6,"color",[0.3,0.3,0.3])
        hold on
        for b = 1:length(aclick.s100)
            plot([1,1]*aclick.s100(b)/1000,[-1,1],"r")
        end
        xlim([0,60])
        ylim([-1,1])
        xticklabels([])
        yticklabels([])
        axis off
    end
    axes('Units','normalized','Position',[0.1 0.47 0.8 0.25]);
    bar(ratbeat(i,:),'FaceColor',[1 .5 0])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    yticks([-3,0,3])
    ylim([-3.4,3.4])
    xticks([1,132])
    axes('Units','normalized','Position',[0.1 0.1 0.8 0.25]);
    bar(humanbeat(i,:),'FaceColor',[0 .5 1])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    yticks([-3,0,3])
    ylim([-3.4,3.4])
    xticks([1,132])
    %saveas(gca, path1 + "summary3/ratbeat2_"+lnames(i)+akinds(apkup) + ".png")
end

%% plot compared figure of human and rat (histogram)

lnames = ["75%","100%","200%","400%"];

for i = 1:length(lnames)
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.4]);
    axes('Units','normalized','Position',[0.1 0.47 0.8 0.25]);
    histogram(ratbeat_(i,:)*9.8,20, "normalization","probability",'FaceColor',[0.8 .4 0])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    if apkup == 3
        xlim([8,40]*10)
        xticks(100:100:500)
    elseif apkup == 1
        xlim([10.5,16])
    end
    ylim([0,0.2])
    axes('Units','normalized','Position',[0.1 0.1 0.8 0.25]);
    histogram(humanbeat_(i,:)*9.8,20, "normalization","probability",'FaceColor',[0 .4 0.8])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    if apkup == 3
        xlim([0,5]*10)
    elseif apkup == 1
        xlim([10.5,12])
    end
    ylim([0,0.2])
    %saveas(gca, path1 + "summary3/ratbeathist2_"+lnames(i)+akinds(apkup) + ".png")
end

%% correlation of human and rat

lnames = ["75%","100%","200%","400%"];

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.11, 1]);
for i = 1:length(lnames)
    subplot(4,1,i)
    mdl = fitlm(ratbeat(i,:), humanbeat(i,:), "Exclude", ratbeat(i,:) > 2.5);
    coef = table2array(mdl.Coefficients);
    disp([lnames(i),coef(2,1),coef(2,end)])
    scatter(ratbeat(i,:),humanbeat(i,:), "k","fill")
    hold on
    if i == 2
        plot([-3,3],[-3,3]*coef(2,1)+coef(1,1),"r", "linewidth",1)
    end
    xlim([-3,3])
    ylim([-3,3])
    yticks([-2,0,2])
    ax = gca;
    ax.YAxis.FontSize = 11;
    ax.XAxis.FontSize = 11;
    crr = corrcoef(ratbeat(i,:), humanbeat(i,:));
end

%saveas(gca, path1 + "summary3/corrbeat_rathuman2_"+akinds(apkup) + ".png")


%% calcurate ratbeat

lnames = ["75%","100%","200%","400%"];
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:4
    subplot(4,1,i)
    bar(ratbeat(i,:))
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    yticks([-3,3])
    ylim([-3,3])
    xticks([1,132])
end

xlabel("Beat Index", "Fontsize",23)
saveas(gca, path1 + "summary3/ratbeat2_"+akinds(apkup) + ".png")


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:length(lnames)
    subplot(4,1,i)
    histogram(ratbeat_(i,:),20, "normalization","probability")
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    xlim([min(min(ratbeat_)),max(max(ratbeat_))])
    ylim([0,0.2])
end

xlabel(akinds(apkup) + "", "Fontsize",23)
%saveas(gca, path1 + "summary3/ratbeathist2_"+akinds(apkup) + ".png")

%% plot rat only

ratfeat_ = squeeze(nanmean(beatfeatr));
ratfeat = zscore(ratfeat_,0,2);
lnames = ["75%","100%","200%","400%"];
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:4
    subplot(4,1,i)
    bar(ratfeat(i,:))
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    yticks([-3,3])
    ylim([-3,3])
    xticks([1,132])
end
xlabel("Beat Index", "Fontsize",23)
%saveas(gca, path1 + "summary3/ratfeat2_"+akinds(apkup) + ".png")


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:length(lnames)
    subplot(4,1,i)
    histogram(ratfeat_(i,:),20)
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    xlim([min(min(ratfeat_)),max(max(ratfeat_))])
end

xlabel(akinds(apkup) + "", "Fontsize",23)
%saveas(gca, path1 + "summary3/ratfeathist2_"+akinds(apkup) + ".png")


%% calcurate humanbeat

lnames = ["75%","100%","200%","400%"];

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:4
    subplot(4,1,i)
    bar(humanbeat(i,:))
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    yticks([-3,3])
    ylim([-3,3])
    xticks([1,132])
end

xlabel("Beat Index", "Fontsize",23)
%saveas(gca, path2 + "summary3/humanbeat2_"+akinds(apkup) + ".png")


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:length(lnames)
    subplot(4,1,i)
    histogram(humanbeat_(i,:),20, "normalization","probability")
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    xlim([min(min(humanbeat_)),max(max(humanbeat_))])
    ylim([0,0.2])
end

xlabel(akinds(apkup) + "", "Fontsize",23)
%saveas(gca, path2 + "summary3/humanbeathist2_"+akinds(apkup) + ".png")

%% plot human only

humanfeat_ = squeeze(nanmean(beatfeat));
humanfeat = zscore(humanfeat_,0,2);

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:4
    subplot(4,1,i)
    bar(humanfeat(i,:))
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    yticks([-3,3])
    ylim([-3,3])
    xticks([1,132])
end
xlabel("Beat Index", "Fontsize",23)
%saveas(gca, path2 + "summary3/humanfeat2_"+akinds(apkup) + ".png")

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.3, 0.7]);
for i = 1:length(lnames)
    subplot(4,1,i)
    histogram(humanfeat_(i,:),20)
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    xlim([min(min(humanfeat_)),max(max(humanfeat_))])
end

xlabel(akinds(apkup) + "", "Fontsize",23)
%saveas(gca, path2 + "summary3/humanfeathist2_"+akinds(apkup) + ".png")


%% correlation of ratfeat and humanfeat

lnames = ["75%","100%","200%","400%"];

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.2, 1]);
for i = 1:length(lnames)
    subplot(4,1,i)
    X = [ones(1,132); ratfeat(i,:)]';
    b = X\humanbeat(i,:)';
    disp([lnames(i),b(2)])
    scatter(ratfeat(i,:),humanfeat(i,:), "k","fill")
    hold on
    plot([-3,3],[-3,3]*b(2)+b(1),"r", "linewidth",1)
    xlim([-3,3])
    ylim([-3,3])
    title(lnames(i), "Fontsize",23)
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel("Human", "Fontsize", 20)
    crr = corrcoef(ratbeat(i,:), humanbeat(i,:));
end

xlabel("Rat", "Fontsize", 20)
%saveas(gca, path1 + "summary3/corrfeat_rathuman2_"+akinds(apkup) + ".png")


%% WAV data with beat lines

[yda,Fs] = audioread(spath + "sonata_" + cnames(2) + ".wav");
%x = (1:length(y))/Fs;
xda = (1:length(yda));

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.25]);
plot(xda/Fs,yda(:,1),"color",[0.3,0.3,0.3])
hold on
for b = 1:length(aclick.s100)
    plot([1,1]*aclick.s100(b)/1000,[-0.6,0.6],"r")
end

xlim([0,60])
ylim([-1,1])
xticklabels([])
yticklabels([])
axis off
%[h,icons] = legend("raw waveform","Fontsize",28);
%icons2 = findobj(icons,'Type','line');
%set(icons2,'linewidth',3);

saveas(gca,path1 + "summary3/detect_peak_raw_wthbeat.png")


