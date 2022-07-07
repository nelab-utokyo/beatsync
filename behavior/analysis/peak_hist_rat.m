%% load music data
cnames = ["s75","s100","s200","s400","click","off"];

load('/Users/TP/research/Tlab/code/sound/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);

spath = "/Users/TP/research/Tlab/sound/material/";

scand = ["75","100","200","400","click"];
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

apkup = 3;  %choose the calculation method below

akinds = ["acc","Z(acc)","jerk","Z(jerk)"];
mkinds = ["s75", "s100", "s200", "s400"];


%% calculate peak acc of each beat (rat)

path1 = "/Users/TP/research/Tlab/behavior/acceleration/";
load(path1 + "dno.mat")
allfiles = [dno1,dno2,dno3,dno4,dno5,dno6,dno7,dno8,dno9,dno11];

for i = 1:size(allfiles,1)
    for j = 1:size(allfiles,2)
        allfiles(i,j).folder = replace(allfiles(i,j).folder,"yoshiki","TP");
    end
end

beatfeatr = zeros(size(allfiles,2),length(mkinds),132);  %collecting the acc value on each beat
maxbeatr = nan(size(allfiles,2),length(mkinds),132,2); %collecting the maximum acc between beats
peakhist = struct();
certify = zeros(size(allfiles,2),length(mkinds));

for k = 1:size(allfiles,2)
    files = allfiles(1,k);
    for j = 1:length(files)
        load(files(j).folder + "/" + files(j).name)
        dt0 = mean(diff(stamp));
        disp(files(j).name)
        disp([dt0,var(diff(stamp))]*1000)
        monoff = ones(size(ontime,1),5);
        for i =1:size(ontime,1)
            [M,Ion] = min(abs(stamp - ontime(i,1) + ontime(i,2)/1000));
            [M,Ioff] = min(abs(stamp - offtime(i,1)- (mlns1(offtime(i,3)+1)-offtime(i,2))/1000));
            monoff(i,2) = Ion;
            monoff(i,3) = Ioff;
            monoff(i,4) = ontime(i,3)+1;
            if i < size(ontime,1)
                monoff(i+1, 1) = Ioff + 1;
            end
            peakhist.("a"+k).(mkinds(monoff(i,4))) = [];
        end
        if apkup == 1
            tmp2_ = acc;
        elseif apkup == 2
            tmp2_ = zscore(acc,0,2);
        elseif apkup == 3
            tmp2_ = diff(acc,1,2)/dt0;
            stamp = (stamp(1:end-1) + stamp(2:end))/2; 
        elseif apkup == 4
            tmp2_ = zscore(diff(acc,1,2)/dt0,0,2);
            stamp = (stamp(1:end-1) + stamp(2:end))/2;
        end
        tmp2 = sqrt(tmp2_(1,:).^2+tmp2_(2,:).^2+tmp2_(3,:).^2);
        %
        for i = 1:size(ontime,1)
            st1 = monoff(i,2);
            Xst = ontime(i,1) - ontime(i,2)/1000;
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
                %tmp = (stamp >= (rng(1))) & (stamp < (rng(2)));
                tmp = (stamp >= (rng(1)-dt/2)) & (stamp <= (rng(1)+dt/2));
                [~,idt] = min(abs(stamp - rng(1)));
                beatfeatr(k,mpkup,rem(i2,132)+1) = beatfeatr(k,mpkup,rem(i2,132)+1) + tmp2(idt)/2;
                if ~isempty(tmp2(tmp))
                    [Mr, Ir] = max(tmp2(tmp));
                    maxbeatr(k,mpkup,i2,1) = maxbeatr(k,mpkup,i2,1) + Mr;
                    maxbeatr(k,mpkup,i2,2) = maxbeatr(k,mpkup,i2,2) + 1;
                    dif = [(stamp(tmp) - rng(1))/dt;tmp2(tmp);tmp2(tmp)*0+i2];
                    peakhist.("a"+k).(mkinds(monoff(i,4))) = [peakhist.("a"+k).(mkinds(monoff(i,4))),dif];      
                    %{
                    fnd = find(tmp);
                    [Mr, Ir] = max(tmp2(fnd(1)-1:fnd(end)+1));
                    if Ir > 1 && Ir < length(fnd)+2
                        maxbeatr(k,mpkup,i2,3+(i>4)) = (stamp(fnd(Ir-1)) - rng(1))/range(rng);
                    else
                        maxbeatr(k,mpkup,i2,3+(i>4)) = NaN;
                    end
                    %}
                end
            end
        end 
        %}
    end
end

maxbeatr(:,:,:,1) = maxbeatr(:,:,:,1)./maxbeatr(:,:,:,2);

thresh = nanmean(maxbeatr(:,:,:,1),3) + 2.5*nanstd(maxbeatr(:,:,:,1),0,3);
maxbeatr(maxbeatr(:,:,:,1) - thresh > 0) = NaN;
ratbeat_ = squeeze(nanmean(maxbeatr(:,:,:,1)));
ratbeat = (ratbeat_ - nanmean(ratbeat_,2))./nanstd(ratbeat_,0,2);


%%

path1 = "/Users/TP/research/Tlab/behavior/acceleration/";

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 1, 0.4]);
ttls = [75,100,200,400] + "%";

nbin = 9;
nbin2 = 7;
for mnum = 1:length(ttls)
    distrbx = [];
    distrby = [];
    for k = 1:size(allfiles,2)
        if sum([1,3,6,8,10] == k) == 0
            continue
        end
        distrbx = [distrbx, peakhist.("a"+k).(mkinds(mnum))(1,:)];
        distrby = [distrby, zscore(peakhist.("a"+k).(mkinds(mnum))(2,:))]; 
    end
    
    distrby((distrbx < -0.5)|(distrbx > 0.5)) = [];
    distrbx((distrbx < -0.5)|(distrbx > 0.5)) = [];
    
    subplot(1, length(ttls), mnum)
    [histw, interval] = histwc(distrbx, distrby, nbin);
    [histwx, interval] = histwc(distrbx, distrby*0+1, nbin);
    histw = histw./histwx;
    distrbx_ = [distrbx-1, distrbx, distrbx+1];
    distrby_ = [distrby, distrby, distrby];
    [f,xi] = kernel(distrbx_, distrby_, nbin2*3);
    ftmp = f;
    ftmp((xi < -0.5)|(xi > 0.5)) = 0;
    [Mf, If] = max(ftmp);
    bar(interval+0.5/nbin, histw)
    hold on
    plot(xi, f/length(distrbx)*nbin2, "Color",[0.7,0.3,0], "linewidth",1.5)
    %plot([1,1]*xi(If), [-1,20],"--","Color",[0.5,0.5,0],"linewidth",2)
    %plot([-1,1], [1,1]*mean(distrby),"g--","linewidth",2)
    interval = linspace(-0.5, 0.5, nbin+1);
    for i = 1:length(interval)-1
        tmp = distrby(distrbx >= interval(i) & distrbx < interval(i+1));
        [h,pv] = ttest(tmp, mean(distrby),"Tail","right");
        if pv < 0.05
            disp([i,pv])
            %plot(interval(i)+0.5/nbin, 0.05, "k*","markersize",20,"linewidth",2)
        end
    end
    xticks(-0.5:0.25:0.5)
    xticklabels(["-π","-π/2","0","π/2","π"])
    xlim([-0.5,0.5])
    ylim([-0.06,0.07])
    yticks(-0.05:0.05:0.5)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    title(ttls(mnum),"Fontsize",25)
    if mnum == 1
        ylabel("Zscore","Fontsize",25)
    end
end
%saveas(gca, path1 + "summary4/peakhist.png")


%%

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 1, 0.4]);
ttls = "pos."+ [1,2,3,4];

nbin = 9;
nbin2 = 7;
mnum = 2;
distrbx = [];
distrby = [];
distrbz = [];
for k = 1:size(allfiles,2)
    if sum([1,3,6,8,10] == k) == 0
        continue
    end
    disp(k)
    distrbx = [distrbx, peakhist.("a"+k).(mkinds(mnum))(1,:)];
    distrby = [distrby, zscore(peakhist.("a"+k).(mkinds(mnum))(2,:))];
    distrbz = [distrbz, peakhist.("a"+k).(mkinds(mnum))(3,:)];
end
distrby((distrbx < -0.5)|(distrbx > 0.5)) = [];
distrbz((distrbx < -0.5)|(distrbx > 0.5)) = [];
distrbx((distrbx < -0.5)|(distrbx > 0.5)) = [];

for p = 1:1 %length(ttls)
    subplot(1, length(ttls), p)
    distx = distrbx(rem(distrbz-1,4)==(p-1));
    disty = distrby(rem(distrbz-1,4)==(p-1));
    [histw, interval] = histwc(distx, disty, nbin);
    [histwx, interval] = histwc(distx, disty*0+1, nbin);
    histw = histw./histwx;
    distrbx_ = [distx-1, distx, distx+1];
    distrby_ = [disty, disty, disty];
    [f,xi] = kernel(distrbx_, distrby_, nbin2*3);
    ftmp = f;
    ftmp((xi < -0.5)|(xi > 0.5)) = 0;
    [Mf, If] = max(ftmp);
    bar(interval+0.5/nbin, histw)
    hold on
    plot(xi, f/length(distx)*nbin2, "Color",[0.7,0.3,0], "linewidth",1.5)
    %plot([1,1]*xi(If), [-1,20],"--","Color",[0.5,0.5,0],"linewidth",2)
    %plot([-1,1], [1,1]*mean(distrby),"g--","linewidth",2)
    xticks(-0.5:0.25:0.5)
    xticklabels(["-π","-π/2","0","π/2","π"])
    xlim([-0.5,0.5])
    ylim([-0.1,0.25])
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    title(ttls(p),"Fontsize",25)
    if p == 1
        ylabel("Zscore","Fontsize",25)
    end
    interval = linspace(-0.5, 0.5, nbin+1);
    for i = 1:length(interval)-1
        tmp = disty(distx >= interval(i) & distx < interval(i+1));
        [h,pv] = ttest(tmp, mean(distrby),"Tail","right");
        if pv < 0.05
            disp([i,pv, mean(tmp), mean(distrby)])
            [p,h] = ranksum(tmp, tmp*0 + mean(distrby),"Tail","both")
            plot(interval(i)+0.5/nbin, 0.2, "k*","markersize",20,"linewidth",2)
        end
    end
    yticks(-0.1:0.1:0.5)
end

%saveas(gca, path1 + "summary4/peakhist_" +scand(mnum) + ".png")



%% calculate peak acc of each beat (human)

path2 = "/Users/TP/research/Tlab/behavior/human/";
load(path2 + "allfiles.mat")
allfiles([2,6,10]) = [];
for i = 1:size(allfiles,1)
    for j = 1:size(allfiles,2)
        allfiles(i,j).folder = replace(allfiles(i,j).folder,"yoshiki","TP");
    end
end

mkinds = ["s75", "s100", "s200", "s400", "click"];
akinds = ["acc","Z(acc)","jerk","Z(jerk)"];

beatfeat = zeros(size(allfiles,2),length(mkinds),132);
maxbeat = zeros(size(allfiles,2),length(mkinds),132,2);
peakhisth = struct();

for k = 1:size(allfiles,2)
    files = allfiles(1,k);
    for j = 1:length(files)
        load(files(j).folder + "/" + files(j).name)
        dt0 = mean(diff(stamp));
        disp(files(j).name)
        disp([dt0,var(diff(stamp))]*1000)
        monoff = ones(size(ontime,1),5);
        for i =1:size(ontime,1)
            [M,Ion] = min(abs(stamp - ontime(i,1)+ ontime(i,2)/1000));
            [M,Ioff] = min(abs(stamp - offtime(i,1)- (mlns1(offtime(i,3)+1)-offtime(i,2))/1000));
            monoff(i,2) = Ion;
            monoff(i,3) = Ioff;
            monoff(i,4) = ontime(i,3)+1;
            if i < size(ontime,1)
                monoff(i+1, 1) = Ioff + 1;
            end
            peakhisth.("a"+k).(mkinds(monoff(i,4))) = [];
        end
        if apkup == 1
            tmp2_ = acc;
        elseif apkup == 2
            tmp2_ = zscore(acc,0,2);
        elseif apkup == 3
            tmp2_ = diff(acc,1,2)/dt0;
            stamp = (stamp(1:end-1) + stamp(2:end))/2;
        elseif apkup == 4
            tmp2_ = zscore(diff(acc,1,2)/dt0,0,2);
            stamp = (stamp(1:end-1) + stamp(2:end))/2;
        end
        tmp2 = sqrt(tmp2_(1,:).^2+tmp2_(2,:).^2+tmp2_(3,:).^2);
        for i = 1:size(ontime,1)
            st1 = monoff(i,2);
            Xst = ontime(i,1) - ontime(i,2)/1000;
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
                tmp = (stamp >= rng(1)-dt/2) & (stamp < rng(1)+dt/2);
                [~,idt] = min(abs(stamp - rng(1)));
                beatfeat(k,mpkup,i2) = beatfeat(k,mpkup,i2) + tmp2(idt)/2;
                if ~isempty(tmp2(tmp))
                    [Mr, Ir] = max(tmp2(tmp));
                    maxbeat(k,mpkup,i2,1) = maxbeat(k,mpkup,i2,1) + Mr;
                    maxbeat(k,mpkup,i2,2) = maxbeat(k,mpkup,i2,2) + 1;
                    dif = [(stamp(tmp) - rng(1))/dt;tmp2(tmp);tmp2(tmp)*0+i2];
                    peakhisth.("a"+k).(mkinds(monoff(i,4))) = [peakhisth.("a"+k).(mkinds(monoff(i,4))),dif];
                    %fnd = find(tmp);
                    %[Mr, Ir] = max(tmp2(fnd(1)-1:fnd(end)+1));
                    %if Ir > 1 && Ir < length(fnd)+2
                    %    maxbeat(k,mpkup,i2,3) = (stamp(fnd(Ir-1)) - rng(1))/range(rng);
                    %else
                    %R    maxbeat(k,mpkup,i2,3) = NaN;
                    %end
                end
            end
        end 
    end
end

maxbeat(:,:,:,1) = maxbeat(:,:,:,1)./maxbeat(:,:,:,2);
humanbeat_ = squeeze(nanmean(maxbeat(:,:,:,1)));
humanbeat = zscore(humanbeat_,0,2);


%%

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 1, 0.4]);
ttls = [75,100,200,400] + "%";

nbin = 9;
nbin2 = 9;
for mnum = 1:length(ttls)
    distrbx = [];
    distrby = [];
    for k = 1:size(allfiles,2)
        disp(k)
        distrbx = [distrbx, peakhisth.("a"+k).(mkinds(mnum))(1,:)];
        distrby = [distrby, zscore(peakhisth.("a"+k).(mkinds(mnum))(2,:))];
    end
    distrby((distrbx < -0.5)|(distrbx > 0.5)) = [];
    subplot(1, length(ttls), mnum)
    [histw, interval] = histwc(distrbx, distrby, nbin);
    [histwx, interval] = histwc(distrbx, distrby*0+1, nbin);
    histw = histw./histwx;
    distrbx_ = [distrbx-1, distrbx, distrbx+1];
    distrby_ = [distrby, distrby, distrby];
    [f,xi] = kernel(distrbx_, distrby_, nbin2*3);
    ftmp = f;
    ftmp((xi < -0.5)|(xi > 0.5)) = 0;
    [Mf, If] = max(ftmp);
    bar(interval+0.5/nbin, histw)
    hold on
    plot(xi, f/length(distrbx)*nbin2, "Color",[0.7,0.3,0], "linewidth",1.5)
    %plot([1,1]*xi(If), [-1,1.5],"--","Color",[0.5,0.5,0],"linewidth",2)
    %plot([-1,1], [1,1]*mean(distrby),"g--","linewidth",2)
    interval = linspace(-0.5, 0.5, nbin+1);
    for i = 1:length(interval)-1
        tmp = distrby(distrbx >= interval(i) & distrbx < interval(i+1));
        [h,pv] = ttest(tmp, mean(distrby),"Tail","right");
        if pv < 0.05
            disp([i,pv])
            %plot(interval(i)+0.5/nbin, 0.4, "k*","markersize",20,"linewidth",2)
        end
    end
    xticks(-0.5:0.25:0.5)
    xticklabels(["-π","-π/2","0","π/2","π"])
    xlim([-0.5,0.5])
    ylim([-0.35,0.5])
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    title(ttls(mnum),"Fontsize",25)
    if mnum == 1
        ylabel("Zscore","Fontsize",25)
    end
    
    yticks(-1:0.2:1)
end
saveas(gca, path1 + "summary4/peakhist_human.png")


%%

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 1, 0.4]);
ttls = "pos."+ [1,2,3,4];

mnum = 2;
distrbx = [];
distrby = [];
distrbz = [];
for k = 1:size(allfiles,2)
    distrbx = [distrbx, peakhisth.("a"+k).(mkinds(mnum))(1,:)];
    distrby = [distrby, zscore(peakhisth.("a"+k).(mkinds(mnum))(2,:))];
    distrbz = [distrbz, peakhisth.("a"+k).(mkinds(mnum))(3,:)];
end
distrby((distrbx < -0.5)|(distrbx > 0.5)) = [];
distrbz((distrbx < -0.5)|(distrbx > 0.5)) = [];
distrbx((distrbx < -0.5)|(distrbx > 0.5)) = [];

for p = 1:length(ttls)
    subplot(1, length(ttls), p)
    distx = distrbx(rem(distrbz-1,4)==(p-1));
    disty = distrby(rem(distrbz-1,4)==(p-1));
    [histw, interval] = histwc(distx, disty, nbin);
    [histwx, interval] = histwc(distx, disty*0+1, nbin);
    histw = histw./histwx;
    distrbx_ = [distx-1, distx, distx+1];
    distrby_ = [disty, disty, disty];
    [f,xi] = kernel(distrbx_, distrby_, nbin2*3);
    ftmp = f;
    ftmp((xi < -0.5)|(xi > 0.5)) = 0;
    [Mf, If] = max(ftmp);
    
    bar(interval+0.5/nbin, histw)
    hold on
    plot(xi, f/length(distx)*nbin2, "Color",[0.7,0.3,0], "linewidth",1.5)
    %plot([1,1]*xi(If), [-1,1.5],"--","Color",[0.5,0.5,0],"linewidth",2)
    %plot([-1,1], [1,1]*mean(distrby),"g--","linewidth",2)
    interval = linspace(-0.5, 0.5, nbin+1);
    for i = 1:length(interval)-1
        tmp = disty(distx >= interval(i) & distx < interval(i+1));
        [h,pv] = ttest(tmp, mean(distrby),"Tail","right");
        if pv < 0.05
            disp([i,pv])
            plot(interval(i)+0.5/nbin, 0.4, "k*","markersize",20,"linewidth",2)
        end
    end
    xticks(-0.5:0.25:0.5)
    xticklabels(["-π","-π/2","0","π/2","π"])
    xlim([-0.5,0.5])
    ylim([-0.35,0.5])
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    title(ttls(p),"Fontsize",25)
    if p == 1
        ylabel("Zscore","Fontsize",25)
    end
    yticks(-1:0.2:1)
end

saveas(gca, path1 + "summary4/peakhist_" +scand(mnum) + "_human.png")




%%
function [histw, interval] = histwc(vv, ww, nbins)
  minV  = -0.5;
  maxV  = 0.5;
  vinterval = linspace(minV, maxV, nbins+1);
  interval = linspace(minV, maxV, nbins+1);
  interval(end) = [];
  histw = zeros(nbins, 1);
  for i=1:length(vv)
    ind = find(vinterval < vv(i), 1, 'last');
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
end

function [f,xi] = kernel(data, weight, nbins)
  minV  = -1.5;
  maxV  = 1.5;
  delta = (maxV-minV)/nbins;
  xi = linspace(minV, maxV, nbins);   
  f = zeros(nbins, 1);
  for i = 1:nbins
      tmpx = (xi(i) - data)/delta;
      k = exp(-(tmpx.^2)/2)/sqrt(2*pi);
      f(i) = dot(k,weight);
  end
end

