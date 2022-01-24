paths = ["20200908","20200910","20200915","20200917","20201006"] + "/";

%% music response

%load(path + "shuf_param.mat")

load("/Users/yoshiki/research/Tlab/code/sound/trigger/click_trg2.mat")
cands = [75,100,150,200,300,400];
start = 1505;

TC = 1:96;

%%
%%%%%%%%%%%%%%%%%%%%%%
%evoked response music
%%%%%%%%%%%%%%%%%%%%%%%

path = paths(1);

load(path + "/processed/exdata_ex2.mat")
resm = struct();

for c = 1:length(cands)
    data = df.("c"+cands(c));
    trg = round(clspeed.("s"+cands(c))*1000 + start);
    [resp,reslfp] = data2matrix3(data,trg,TC,0);
    resp = squeeze(mean(resp,2));
    resm.("c"+cands(c)) = resp;
    reslfp = squeeze(mean(reslfp,2));
    resm.("cl"+cands(c)) = reslfp;
end
save(path + "buildmat/evoke_resm3.mat","resm")

%% spontaneous evoked response
%load(path + "shuf_param.mat")
TC = 1:96;
load(path + "processed/exdata_ex2_spon.mat")

lng = 10900;
df.trg = 1:lng;
resp = zeros(length(TC),10,length(df.trg)+1);
reslfp = zeros(length(TC),10,length(df.trg));

for e = 1:length(TC)
    elenum = TC(e);
    for i = 0:9
        temp = df.("s"+i){elenum};
        temp(temp > length(df.trg)) = [];
        resp(e,i+1,temp+1) = resp(e,i+1,temp+1) + 1;
        reslfp(e,i+1,:) = df.("d"+i)(elenum,1:lng);
    end
end

n = 1000;
temp = zeros(length(TC),n);
temp2 = zeros(length(TC),n);

resp_ = squeeze(mean(resp, 2))*1000;
xs = round(rand(n,1)*(lng-50));
for i = 1:n
    temp(:,i) = mean(resp_(:,(xs(i)+25):(xs(i)+50)),2) - mean(resp_(:,(xs(i)+15):(xs(i)+24)),2);
    temp2(:,i) = mean(resp_(:,(xs(i)+25):(xs(i)+50)),2);
end

spntevk = [mean(temp,2),std(temp,0,2)];
spntevk2 = [mean(temp2,2),std(temp2,0,2)];

save(path + "buildmat/spnt_evok.mat","spntevk","spntevk2")


%%
%%%%%%%%%%%%%%%%%%
%evoked response p
%%%%%%%%%%%%%%%%%%

%{
%load(path + "shuf_param.mat")
TC = 1:96;
load(path + "processed/exdata_p2p.mat")
names = "p"+[1,2,4,8,16,32];
namesl = "pl"+[1,2,4,8,16,32];
evp = struct();

for n = 1:length(names)
    name = names(n) + "_15";
    data = df.(name);
    [resp,reslfp] = data2matrix2(data,data.trg2,TC,0);
    resp = squeeze(mean(resp,2));
    evp.(names(n)) = resp;
    reslfp = squeeze(mean(reslfp,2));
    evp.(namesl(n)) = reslfp;
end

save(path + "buildmat/evoke_res2.mat","evp")

evp = struct();

for n = 1:length(names)
    name = names(n) + "_15";
    data = df.(name);
    [resp,reslfp] = data2matrix3(data,data.trg2,TC,0);
    resp = squeeze(mean(resp,2));
    evp.(names(n)) = resp;
    reslfp = squeeze(mean(reslfp,2));
    evp.(namesl(n)) = reslfp;
end

save(path + "buildmat/evoke_res3.mat","evp")



%%
%%%%%%%%%%%%%%%%%%
%evoked response f
%%%%%%%%%%%%%%%%%%

%load(path + "shuf_param.mat")
TC = 1:96;
load(path + "processed/exdata_p2f.mat")
names = "f"+[4,8,16,32];
namesl = "fl"+[4,8,16,32];
evf = struct();

for n = 1:length(names)
    name = names(n) + "_15";
    data = df.(name);
    [resp,reslfp] = data2matrix3(data,data.trg2,TC,0);
    resp = squeeze(mean(resp,2));
    evf.(names(n)) = resp;
    reslfp = squeeze(mean(reslfp,2));
    evf.(namesl(n)) = reslfp;
end

save(path + "buildmat/evoke_res3.mat","evf","-append")

%}