paths = ["20200908","20200910","20200915","20200917","20201006","20210106","20210107"] + "/";
%paths = ["20210106","20210107"] + "/";

%% music response

%load(path + "shuf_param.mat")

load('/Users/TP/research/Tlab/code/sound/click_trg2.mat')
cands = [75,100,150,200,300,400];
start = 1505;

TC = 1:96;

for path = paths
    %
    load(path + "/processed/exdata_ex1.mat")
    disp(path)
    resm = struct();
    
    for c = 1:length(cands)
        data = df.("c"+cands(c));
        trg = round(clspeed.("s"+cands(c))*1000 + start);
        [resp,reslfp,ssum] = data2matrix3(data,trg,TC,0);
        %[resp,reslfp] = data2matrix(data,trg,TC,0);
        resp = squeeze(mean(resp,2));
        resm.("c"+cands(c)) = resp;
        %resm.("s"+cands(c)) = ssum;
        reslfp = squeeze(mean(reslfp,2));
        resm.("cl"+cands(c)) = reslfp;
    end
    
    save(path + "buildmat/evoke_resm3.mat","resm")
    %}
    try
        load(path + "processed/exdata_p2_spon.mat")
    catch
        load(path + "processed/exdata_ex2_spon.mat")
    end
    lng = 10900;
    df.trg = 1:lng;
    resp = zeros(length(TC),10,length(df.trg)+1);
    reslfp = zeros(length(TC),10,length(df.trg));
    
    for e = 1:length(TC)
        elenum = TC(e);
        for i = 0:9
            try
                reslfp(e,i+1,:) = df.("d"+i)(elenum,1:lng);
                temp = df.("s"+i){elenum};
            catch
                reslfp(e,i+1,:) = df.csp.("d"+(i+1))(elenum,1:lng);
                temp = df.csp.("s"+(i+1)){elenum};
            end
            temp(temp > length(df.trg)) = [];
            resp(e,i+1,temp+1) = resp(e,i+1,temp+1) + 1;
        end
    end
    
    n = 1000;
    temp = zeros(length(TC),n);
    temp2 = zeros(length(TC),n);
    
    resp_ = squeeze(mean(resp, 2))*1000;
    xs = round(rand(n,1)*(lng-50));
    for i = 1:n
        temp(:,i) = mean(resp_(:,(xs(i)+25):(xs(i)+50)),2) - mean(resp_(:,(xs(i)+20):(xs(i)+24)),2);
        temp2(:,i) = mean(resp_(:,(xs(i)+25):(xs(i)+50)),2);
    end
    
    spntevk = [mean(temp,2),std(temp,0,2)];
    spntevk2 = [mean(temp2,2),std(temp2,0,2)];
    
    save(path + "buildmat/spnt_evok.mat","spntevk","spntevk2")

end



%% individual data processing
%
%%%%%%%%%%%%%%%%%%
%evoked response 
%%%%%%%%%%%%%%%%%%

path = "20200902/";
load(path + "/processed/exdata_ex1.mat")
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
%}

%% spontaneous evoked response

path = "20200902/";
load(path + "shuf_param.mat")
TC = 1:96;
load(path + "processed/exdata_p2_spon.mat")

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
    temp(:,i) = mean(resp_(:,(xs(i)+25):(xs(i)+50)),2) - mean(resp_(:,(xs(i)+20):(xs(i)+24)),2);
    temp2(:,i) = mean(resp_(:,(xs(i)+25):(xs(i)+50)),2);
end

spntevk = [mean(temp,2),std(temp,0,2)];
spntevk2 = [mean(temp2,2),std(temp2,0,2)];

save(path + "buildmat/spnt_evok.mat","spntevk","spntevk2")
%}

%%
%%%%%%%%%%%%%%%%%%
%evoked response p
%%%%%%%%%%%%%%%%%%

path = "20200917/"
load(path + "shuf_param.mat")
TC = 1:96;
load(path + "processed/exdata_p2p.mat")
names = "p"+[1,2,4,8,16,32];
namesl = "pl"+[1,2,4,8,16,32];
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

if exist(path + "buildmat/evoke_res3.mat", 'file')
    save(path + "buildmat/evoke_res3.mat","evp","-append")
else
    save(path + "buildmat/evoke_res3.mat","evp")
end

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

if exist(path + "buildmat/evoke_res2.mat", 'file')
    save(path + "buildmat/evoke_res2.mat","evp","-append")
else
    save(path + "buildmat/evoke_res2.mat","evp")
end

%}
%%
%%%%%%%%%%%%%%%%%%
%evoked response f
%%%%%%%%%%%%%%%%%%

path = "20200820/";
load(path + "shuf_param.mat")
TC = 1:96;
load(path + "processed/exdata_p2f.mat")
names = "f"+[4,8,16,32];
namesl = "fl"+[4,8,16,32];
evf = struct();
evf2 = struct();

for n = 1:length(names)
    name = names(n) + "_15";
    data = df.(name);
    [resp,reslfp] = data2matrix3(data,data.trg2,TC,0);
    resp = squeeze(mean(resp,2));
    evf.(names(n)) = resp;
    reslfp = squeeze(mean(reslfp,2));
    evf.(namesl(n)) = reslfp;
    trg2 = data.trg2;
    trg2 = [trg2, trg2(3:3:end) + round(mean(trg2(2:3:end)- trg2(1:3:end)))]; %include pos. 4
    trg2 = sort(trg2);
    [resp,reslfp] = data2matrix3(data,trg2,TC,0);
    resp = squeeze(mean(resp,2));
    evf2.(names(n)) = resp;
end

if exist(path + "buildmat/evoke_res3.mat", 'file')
    save(path + "buildmat/evoke_res3.mat","evf","evf2","-append")
else
    save(path + "buildmat/evoke_res3.mat","evf","evf2")
end


evf = struct();

for n = 1:length(names)
    name = names(n) + "_15";
    data = df.(name);
    [resp,reslfp] = data2matrix2(data,data.trg2,TC,0);
    resp = squeeze(mean(resp,2));
    evf.(names(n)) = resp;
    reslfp = squeeze(mean(reslfp,2));
    evf.(namesl(n)) = reslfp;
end

if exist(path + "buildmat/evoke_res2.mat", 'file')
    save(path + "buildmat/evoke_res2.mat","evf","-append")
else
    save(path + "buildmat/evoke_res2.mat","evf")
end
%