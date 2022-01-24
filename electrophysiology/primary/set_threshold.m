%%
%%%%%%%%%%%%%%%%%%
%failer response of p
%%%%%%%%%%%%%%%%%%
paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006"];

st = 1;
names = "p"+[1,2,4,8,16,32];

for p = 1:length(paths)
    load("./" + paths(p) + "/buildmat/evoke_res2.mat")
    load("./" + paths(p) +"/buildmat/spnt_evok.mat")
    load("./" + paths(p) + "/shuf_param.mat")
    failrs = zeros(length(TC), length(names), 2);
    rmthresh = (spntevk2(TC,1)+ spntevk2(TC,2)*2);
    failthresh = (spntevk2(TC,1)+ spntevk2(TC,2));
    spn_ = (spntevk2(TC,1)+ spntevk2(TC,2))*0;
    
    for n = 1:length(names)
        name = names(n);
        data = evp.(name)(TC,:);
        data_ = evp.(name)(TC,:);
        if n == 1
            rmv = find(data(:,1) < rmthresh); 
        end
        failrs(:,n,1) = mean(data_ < failthresh,2);
        data_(data_ < failthresh) = NaN;
        failrs(:,n,2) = nanmean(data_,2);
    end
    disp(length(rmv))
    disp(rmv.')
    save("./" + paths(p) + "/buildmat/responsep.mat","rmv","failrs")
end


%%
%%%%%%%%%%%%%%%%%%
%evoked response p
%%%%%%%%%%%%%%%%%%
paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006"];

st = 1;
names = "p"+[1,2,4,8,16,32];

for p = 1:length(paths)
    load("./" + paths(p) + "/buildmat/evoke_res2.mat")
    load("./" + paths(p) +"/buildmat/spnt_evok.mat")
    load("./" + paths(p) + "/shuf_param.mat")
    diffs = zeros(length(TC),length(names),2);
    diffs2 = zeros(length(TC),length(names),2);
    spn_ = (spntevk2(TC,1)+ spntevk2(TC,2))*0;
    
    for n = 1:length(names)
        name = names(n);
        data = evp.(name)(TC,:);
        data_ = evp.(name)(TC,:);
        diffs2(:,n,1) = data(:,1);
        diffs2(:,n,2) = nanmean(data(:,st:end),2);
        data = data - spn_;
        data(data < 0) = 0;
        diffs(:,n,1) = data(:,1);
        diffs(:,n,2) = nanmean(data(:,st:end),2);
    end
    save("./" + paths(p) + "/buildmat/responsep.mat","diffs","diffs2","-append")
end



%%
%%%%%%%%%%%%%%%%%%
%evoked response f
%%%%%%%%%%%%%%%%%%
names = "f"+[4,8,16,32];

for p = 1:length(paths)
    load("./" + paths(p) + "/buildmat/evoke_res3.mat")
    load("./" + paths(p) +"/buildmat/spnt_evok.mat")
    load("./" + paths(p) + "/shuf_param.mat")
    ffs = zeros(length(TC),length(names),2);
    ffs2 = zeros(length(TC),length(names),2);
    spn_ = (spntevk2(TC,1)+ spntevk2(TC,2)*1)*0;
    for n = 1:length(names)
        name = names(n);
        data = evf.(name)(TC,:);
        ffs2(:,n,1) = mean(data(:,1:3:end),2);
        ffs2(:,n,2) = (mean(data(:,2:3:end),2)+(mean(data(:,3:3:end),2)))/2;
        data = data - spn_;
        data(data < 0) = 0;
        st = 1;
        ffs(:,n,1) = mean(data(:,1:3:end),2);
        ffs(:,n,2) = (mean(data(:,2:3:end),2)+(mean(data(:,3:3:end),2)))/2;
    end
    save("./" + paths(p) + "/buildmat/responsef.mat","ffs","ffs2")
end


%%
%%%%%%%%%%%%%%%%%%
%evoked response music
%%%%%%%%%%%%%%%%%%

paths = ["20200908","20200910","20200915","20200917","20201006","20210106","20210107"];
%paths = ["20200908","20200910","20200915","20200917","20201006","20201220","20201222","20210106","20210107","20210108_1","20210108_2"];
cands = [75,100,200,300,400];


load("/Users/yoshiki/research/Tlab/code/sound/trigger/click_trg2.mat")
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];

trg2 = 1:489;
%trg2((diff(clspeed.("s400"))*1000)<25) = NaN;

for p = 1:length(paths)
    load("./"+paths(p)+"/buildmat/evoke_resm3.mat")
    load("./" + paths(p) +"/buildmat/spnt_evok.mat")
    load("./" + paths(p) + "/shuf_param.mat")
    load("./" + paths(p) + "/CF_latency.mat")
    rmv = find(latency(TC) > 20);
    disp(length(TC))
    bres = zeros(length(cands),length(TC),5);
    nbres = zeros(length(cands),length(TC),5);
    spn_ = (spntevk(TC,1)+ spntevk(TC,2))*0;
    for c = 1:length(cands)
        resp_ = resm.("c"+cands(c))(TC,:);
        %resp_ = resp_ - spn_;
        resp_(resp_ < spn_) = 0;
        for i = 1:4
            beats = zeros(1,489);
            beat = onbeat(:,i);
            beat(isnan(beat)) = [];
            beat = trg2(beat);
            beat(isnan(beat)) = []; 
            offbeat = 1:489;
            offbeat(beat) = [];
            offbeat = trg2(offbeat);
            offbeat(isnan(offbeat)) = []; 
            bres(c,:,i) = mean(resp_(:, beat),2);
            nbres(c,:,i) = mean(resp_(:, offbeat),2);
        end
        allbeat = trg2(allbeats);
        allbeat(isnan(allbeat)) = []; 
        offtrg = 1:489;
        offtrg(allbeats) = [];
        offtrg = trg2(offtrg);
        offtrg(isnan(offtrg)) = []; 
        bres(c,:,5) = mean(resp_(:, allbeat),2);
        nbres(c,:,5) = mean(resp_(:, offtrg),2);
    end
    save("./" + paths(p) + "/buildmat/responsem.mat","bres","nbres","rmv")
end


%%
%%%%%%%%%%%%
%create mean response of evp, evf of experimental_data.mat
%%%%%%%%%%%%

paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006"];

plabel = "p" + [1,2,4,8,16,32];
flabel = "f" + [4,8,16,32];

evpsum = struct("p1",[],"p2",[],"p4",[],"p8",[],"p16",[],"p32",[]);
evfsum = struct("f4",[],"f8",[],"f16",[],"f32",[]);
rgns = [];
initres = [];
lastres = [];
lastres2 = [];

for p = 1:length(paths)
    load("./"+paths(p)+"/buildmat/evoke_res2.mat")
    load("./" + paths(p) +"/buildmat/spnt_evok.mat")
    load("./" + paths(p) +"/buildmat/responsep.mat")
    load("./"+paths(p)+"/CF_latency.mat")
    load("./"+paths(p)+"/shuf_param.mat")
    spn_ = (spntevk2(TC,1)+ spntevk2(TC,2))*0;
    spn_(rmv) = [];
    diffs(rmv,:,:) = [];
    tmpa = a1aaf(TC);
    tmpa(rmv) = [];
    rgns = [rgns;tmpa];
    stable_ = [];
    for x = 1:length(plabel)
        tmp = evp.(plabel(x))(TC,:);
        tmp(rmv,:) = [];
        tmp = tmp - spn_;
        tmp(tmp < 0) = 0;
        evpsum.(plabel(x)) = [evpsum.(plabel(x)); tmp];
    end
    initres = [initres; mean(diffs(:,:,1),2)];
    for y = 1:length(flabel)
        tmp = evf.(flabel(y))(TC,:);
        tmp(rmv,:) = [];
        tmp = tmp - spn_;
        tmp(tmp < 0) = 0;
        evfsum.(flabel(y)) = [evfsum.(flabel(y)); tmp]; 
    end
    lastres = [lastres; (mean(tmp(:,(size(tmp,2)*2/3):3:end)+tmp(:,(size(tmp,2)*2/3):3:end),2))/2]; 
    lastres2 = [lastres2; mean(evp.(plabel(end))(TC,320:480),2)];
end


init = [];
for x = 1:length(plabel)
    evpsum.(plabel(x)) = squeeze(nanmean(evpsum.(plabel(x))(find(rgns~=-2),:),1));
    init = [init, evpsum.(plabel(x))(1)];
end
for y = 1:length(flabel)
    evfsum.(flabel(y)) = squeeze(nanmean(evfsum.(flabel(y))(find(rgns~=-2),:),1));
    tmp = mean(evfsum.(flabel(y))(2:3:end)) + mean(evfsum.(flabel(y))(3:3:end));
    disp(tmp/2)
    init = [init, evfsum.(flabel(y))(1)];
end

%save("experiment_data.mat","evpsum","evfsum","init","-append")

%% printing the initial / stable value

lastratio = [];
lastratio2 = [];
for i = 1:5
    if i == 1
        temp = 1:length(lastres);
    elseif i == 2
        temp = find(rgns == 1);
    elseif i == 3
        temp = find(rgns == -1);
    elseif i == 4
        temp = find(rgns ~= 0);
    elseif i == 5
        temp = find(rgns == 0);
    end
    ratio = mean(lastres(temp))/median(initres(temp));
    lastratio = [lastratio, ratio];
    lastratio2 = [lastratio2, mean(lastres(temp))];
end

finalv = 1-lastratio;
disp(finalv)
finalv2 = lastratio2;
disp(finalv2)
%save("experiment_data.mat","finalv","finalv2","-append")


%%

paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006"];

labels = [60,120,240,480];

resf = [];
resp = [];
neunum = [0];
spns = [];
rgns = [];
for p = 1:length(paths)
    load("./"+paths(p)+"/CF_latency.mat")
    load("./"+paths(p)+"/shuf_param.mat")
    load("./"+paths(p)+"/buildmat/responsep2.mat")
    load("./"+paths(p)+"/buildmat/responsef2.mat")
    TC(rmv) = [];
    tmp = a1aaf(TC);
    ffs(rmv,:,:) = [];
    diffs(rmv,:,:) = [];
    resf = [resf; ffs];
    resp = [resp; diffs];
    rgns = [rgns;tmp];
    neunum = [neunum, size(resf,1)];
    disp(size(resp,1))
end

pdata = [mean(resp(:,:,2));std(resp(:,:,2),0,1)/sqrt(size(resp,1))].';
figure()
errorbar(1:size(pdata,1),pdata(:,1),pdata(:,2))

adjusted = zeros(length(labels),2);
adjusted(:,1) = mean(resf(:,:,1),1);
adjusted(:,2) = mean(resf(:,:,2),1);
eradjusted = squeeze(std(resf,0,1)/sqrt(size(resp,1)));
fdata = [adjusted(:,1),eradjusted(:,1),adjusted(:,2),eradjusted(:,2)];

figure()
errorbar(1:length(labels),fdata(:,1),fdata(:,2))
hold on
errorbar(1:length(labels),fdata(:,3),fdata(:,4))

fratio = zeros(length(labels),2);
fratio(:,1) = fdata(:,1)./fdata(:,3);
fratio(:,2) = sqrt((fdata(:,2)./fdata(:,3)).^2 + (fdata(:,4).*fdata(:,1)./fdata(:,3)./fdata(:,3)).^2);

fcontrast_ = zeros(9,4);
for i = 1:length(neunum)-1
    onf = resf(neunum(i)+1:neunum(i+1),:,1);
    offf = resf(neunum(i)+1:neunum(i+1),:,2);
    fcontrast_(i,:) = mean((onf - offf)./(onf + offf));
end
fcontrast = [mean(fcontrast_); std(fcontrast_)/sqrt(9)]';

%fcontrast = zeros(length(labels),2);
%fcontrast(:,1) = mean((resf(:,:,1) - resf(:,:,2))./(resf(:,:,1) + resf(:,:,2)));
%fcontrast(:,2) = std((resf(:,:,1) - resf(:,:,2))./(resf(:,:,1) + resf(:,:,2)))/sqrt(size(resf,1));

figure()
errorbar(1:length(labels),fratio(:,1),fratio(:,2))

figure()
errorbar(1:length(labels),fcontrast(:,1),fcontrast(:,2))

%%

rpdata = zeros(4,6,2);
rfdata = zeros(4,length(labels),4);
rfratio = zeros(4,length(labels),2);
rfcontrast = zeros(4,length(labels),2);

for i = 1:4
    if i == 1
        temp = find(rgns == 1);
    elseif i == 2
        temp = find(rgns == -1);
    elseif i == 3
        temp = find(rgns ~= 0);
    elseif i == 4
        temp = find(rgns == 0);
    end
    rpdata(i,:,:) = [mean(resp(temp,:,2));std(resp(temp,:,2),0,1)/sqrt(length(temp))].';
    adjusted = zeros(length(labels),2);
    adjusted(:,1) = mean(resf(temp,:,1),1);
    adjusted(:,2) = mean(resf(temp,:,2),1);
    eradjusted = squeeze(std(resf(temp,:,:),0,1)/sqrt(length(temp)));
    rfdata_ = [adjusted(:,1),eradjusted(:,1),adjusted(:,2),eradjusted(:,2)];
    rfdata(i,:,:) = rfdata_;
    
    rfratio(i,:,1) = rfdata_(:,1)./rfdata_(:,3);
    rfratio(i,:,2) = sqrt((rfdata_(:,2)./rfdata_(:,3)).^2 + (rfdata_(:,4).*rfdata_(:,1)./rfdata_(:,3)./rfdata_(:,3)).^2);
    
    rfcontrast_ = zeros(9,4);
    for j = 1:length(neunum)-1
        temp_ = temp;
        temp_(temp_ <= neunum(j)) = [];
        temp_(temp_ > neunum(j+1)) = [];
        onf = resf(temp_,:,1);
        offf = resf(temp_,:,2);
        rfcontrast_(j,:) = mean((onf - offf)./(onf + offf));
    end
    rfcontrast(i,:,:) = [mean(rfcontrast_); std(rfcontrast_)/sqrt(9)]';
end

%%
%save("experiment_data.mat","pdata","fdata","fratio","fcontrast","-append")
%save("experiment_data.mat","rpdata","rfdata","rfratio","rfcontrast","-append")

