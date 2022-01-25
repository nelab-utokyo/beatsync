%% load music data
cnames = ["s75","s100","s200","s400","off"];

load('../../electrophysiology/sound/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);

spath = "../software/sound/speed/";

rpnum = [1,1,2,4];
rpnum1 = [1,1,1,4];
ratio = [1,0.9946,1,1];
mlns = [];

aclick = struct();

for i = 1:4
    [y,Fs] = audioread(spath + "sonata_" + cnames(i) + ".wav");
    sound = y(:,1);
    mlns = [mlns, round(length(sound)/Fs*1000)];
    tmp = clspeed.(cnames(i))*1000*ratio(i);
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

%%

cands = ["20201207"] %,"20201208","20201209","20201210","20201211"];
allfiles = [];
for cand = cands
    files = dir(cand+"/*.mat");
    allfile = [];
    for i = 1:1 %length(files)
        if endsWith(files(i).name,".mat")
            allfile = [allfile, files(i)];
        end
    end
    allfiles = [allfiles; allfile];
end

%%
path1 = "/Users/yoshiki/research/Tlab/behavior/acceleration/";
load(path1 + "dno.mat")
allfiles = [dno1,dno2,dno3,dno4,dno5,dno6,dno7,dno8,dno9,dno11];

%%
mkinds = ["s75", "s100", "s200", "s400"];

indvdata = struct();

gap = 0;
onratio = 0.4;
rep = 250;

for k = 1:size(allfiles,2)
    files = allfiles(1:2,k);
    prep = zeros(length(gap),4,2);
    prep2 = zeros(rep,4,2);
    monamp = struct("s75",prep, "s100",prep, "s200",prep, "s400", prep);
    monamp2 = struct("s75",prep2, "s100",prep2, "s200",prep2, "s400", prep2);

    for j = 1:1
        load(files(j).folder + "/" + files(j).name)
        dt = mean(diff(stamp));
        disp(files(j).name)
        disp(dt)
        monoff = ones(size(ontime,1),4);
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
                st = monoff(n,2);
                if i < 4
                    tmp2 = abs(diff(acc(i,st-1:monoff(n,3))));
                else
                    tmp2_ = diff(acc(:,st-1:monoff(n,3)),1,2);
                    tmp2 = sqrt(tmp2_(1,:).^2+tmp2_(2,:).^2+tmp2_(3,:).^2);
                end
                for q = 2:rep+1   %1: correct beat, 2: randomized beat(not used)
                    beatime = round(aclick.(mkinds(monoff(n,4)))/mlns2(monoff(n,4))*length(monoff(n,2):monoff(n,3)));
                    interval = median(diff(beatime(1:132)));                   
                    part = round(onratio*interval);
                    part = part + 1*(part == 0);
                    if q > 1
                        %beatime = beatime + round(rand(size(beatime))*interval);
                        %beatime = round(sort(rand(1,length(beatime))*(beatime(end)-beatime(1))));
                        beatime = randperm(length(tmp2));
                        onb = mean(tmp2(beatime(1:round(onratio*length(beatime)))));
                        offb = mean(tmp2(beatime(round(onratio*length(beatime))+1:end)));
                        monamp2.(mkinds(monoff(n,4)))(q-1,i,:) = monamp2.(mkinds(monoff(n,4)))(q-1,i,:) + reshape([onb,offb],[1,1,2])/2;
                    else
                        g = 1;
                        tmp = zeros(1, length(st:monoff(n,3)));
                        for r = 0:(part-1)
                            r_ = floor(r - (part-3)/2);
                            beatime_ = beatime+gap(g) + r_;
                            beatime_(beatime_<1) = [];
                            tmp(1,beatime_) = 1;
                        end
                        %tmp2(tmp2 > 2) = 2;
                        onb = mean(tmp2(find(tmp==1)));
                        offb = mean(tmp2(find(tmp==0)));
                        monamp.(mkinds(monoff(n,4)))(1,i,:) = monamp.(mkinds(monoff(n,4)))(g,i,:) + reshape([onb,offb],[1,1,2])/2/length(files); 
                    end
                end
            end
        end
    end
    indvdata.("no"+k) = monamp;
    indvdata.("no2"+k) = monamp2;
end


%% calculating the border line

count = 0;
h = [];
dim = 4;

chcker = zeros(length(gap),size(allfiles,2),length(files)*2);

border95 = [];
border99 = [];
for k = 1:size(allfiles,2)
    store = [];
    for knd = 1:4
        tmp = squeeze(indvdata.("no2"+k).(mkinds(knd))(:,dim,:));
        tmp2 = (tmp(:,1) - tmp(:,2))./(tmp(:,1) + tmp(:,2));
        store = [store; tmp2];
    end
    store = sort(store);
    border95 = [border95; store(round(length(store)*0.05)), store(round(length(store)*0.95))];
    border99 = [border99; store(round(length(store)*0.01)), store(round(length(store)*0.99))];
end


%%

save("null_border.mat","border95","border99")
