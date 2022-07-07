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


%% calculate peak acc of each beat (human)

path2 = "/Users/TP/research/Tlab/behavior/human/";
load(path2 + "allfiles_video.mat")
allfiles([2,6,10]) = [];

apkup = 3;
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

subdata = [];
mnum = 2;
nbin = 9;
for k = 1:size(allfiles,2)
    distrbx = peakhisth.("a"+k).(mkinds(mnum))(1,:);
    distrby = zscore(peakhisth.("a"+k).(mkinds(mnum))(2,:));
    distrby(distrbx < -0.5) = [];
    distrby(distrbx > 0.5) = [];
    distrbx(distrbx < -0.5) = [];
    distrbx(distrbx > 0.5) = [];
    [histw, interval] = histwc(distrbx, distrby, nbin);
    [histwx, interval] = histwc(distrbx, distrby*0+1, nbin);
    histw = histw./histwx;
    subdata = [subdata, histw];
end

%%

figure()
newcolors = {'#E00','#E80','#EE0', '#0B0','#0EE','#00E','#50E','#A0E'};
colororder(newcolors)
%interval = linspace(-0.5, 0.5, nbin);

bar(interval+0.5/nbin, subdata', "stack")

xticks(-0.5:0.25:0.5)
xticklabels(["-π","-π/2","0","π/2","π"])
xlim([-0.5,0.5])
yticks(-1:1:1)

ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;
ylabel("Zscore","Fontsize",25)
xlabel('Phase',"Fontsize", 25);

%saveas(gca, "./summary4/jerk_hist_human.png")

%%


xalpha = reshape(repmat((interval+0.5/9)*pi*2,size(allfiles,2),1),1,[]);
yw = reshape(subdata',1,[])
xalpha = xalpha + pi*(yw < 0);
p = circ_rtest(xalpha, abs(yw))  %0.0347


xalpha = reshape(repmat((interval+0.5/9)*pi*2,size(allfiles,2),1),1,[]);
yw = reshape(subdata'/8,1,[]);
xalpha = xalpha + pi*(yw < 0);
p = circ_rtest(xalpha, abs(yw))  %0.7545



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


