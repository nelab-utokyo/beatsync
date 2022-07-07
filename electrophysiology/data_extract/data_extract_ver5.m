%%
%detect start and goal
flag = true;
count = 0;
start = [];
goal = [0];
count_ = 4000;
for i = 1:length(TRIGwave)
    temp = abs(TRIGwave(1,i));
    if flag && temp < 50
        count = count + 1;
        if count == count_ && (i - goal(end)) > 13000 
            goal = [goal,i+500];
        end
    elseif flag && temp > 50 
        if max(TRIGwave(1,i:i+2000)) > 1500
            count_ = 100;
        else
            count_ = 2500;
        end
        flag = false;
        if i > goal(end)+5 && (count < 20000) && length(start) < length(goal)
            start = [start,i-2];
        end
        count = 0;
    elseif ~flag && temp < 50
            flag = true;
    end
end

disp(length(start))  %best value = 60
disp(length(goal))   %the initial 6 values should be correct

start(1,1:length(start)) = start(1,1:length(start)) - 1500;
%%
%start(:,6:6:56) = [];
plot(TRIGwave)
hold on 
plot(start+1500,zeros(length(start)),"ro")
plot(goal,zeros(length(goal),1),"ko")
%plot(goal,zeros(length(goal),1),"ko")
%%

load("shuf_param.mat")
number = 6
remove = [];
dnum_ = [6,6,6,6,6,8];

dnum = dnum_(number);

switch number
    case 1
    shforder = fixed_music1f;
    case 2
    shforder = fixed_music1r;
    case 3
    shforder = fixed_music2f;
    case 4
    shforder = fixed_music2r;
    case 5
    shforder = fixed_music3;
    case 6
    shforder = fixed_music_ex2;
    %name = ["csp","c100","c150","c75","cre","c200","c300","c400"];
    name = ["csp","c100","c150","c75","c50","c200","c300","c400"];
end

if ~isempty(remove)
    shforder(remove) = [];
end

if length(start) ~= length(shforder)
    disp("length does not match")
    pause
end
%%
df = struct();
df.date = "20210108";

for i = 1:dnum
    df.(name(shforder(i)+1)).trg = TRIGwave(start(i):goal(i+1));
end
trigger = [];

for i = 1:dnum
    trigger = [trigger,length(df.(name(i)).trg)];
end
if length(trigger) ~= dnum
    disp("length does not match")
    pause
end

disp(trigger)
%%
counter = zeros(1,dnum);
for i = 1:length(start)
    sss = cell(length(SpikeTiming),1);
    for j = 1:length(SpikeTiming)
        sptm = int32(SpikeTiming{j}/30);
        sptm_ = sptm - start(i);
        sptm_ = sptm_(sptm_ >= 0);
        sptm_ = sptm_(sptm_ < trigger(shforder(i)+1));
        sss{j} = sptm_;
    end
    counter(shforder(i)+1) = counter(shforder(i)+1) + 1;
    df.(name(shforder(i)+1)).("s"+string(counter(shforder(i)+1))) = sss;
    df.(name(shforder(i)+1)).("d"+string(counter(shforder(i)+1))) = analogWaves(:,start(i):start(i)+trigger(shforder(i)+1)-1);
end
disp(counter)
%%

path2_ = "./processed/ex";
if number == 1
    df = spon;
    save(path2_+"data_1fsp.mat","df")
    df = df1;
    save(path2_+"data_1f0.mat","df")
    df = df3;
    save(path2_+"data_1f10.mat","df")
    df = df5;
    save(path2_+"data_1f20.mat","df")
    df = df2;
    save(path2_+"data_1f50.mat","df")
    df = df4;
    save(path2_+"data_1f100.mat","df")
elseif number == 2
    df = spon;
    save(path2_+"data_1rsp.mat","df")
    df = df1;
    save(path2_+"data_1r0.mat","df")
    df = df3;
    save(path2_+"data_1r10.mat","df")
    df = df5;
    save(path2_+"data_1r20.mat","df")
    df = df2;
    save(path2_+"data_1r50.mat","df")
    df = df4;
    save(path2_+"data_1r100.mat","df")
elseif number == 3
    df = spon;
    save(path2_+"data_2fsp.mat","df")
    df = df1;
    save(path2_+"data_2f0.mat","df")
    df = df3;
    save(path2_+"data_2f100p.mat","df")
    df = df5;
    save(path2_+"data_2f100r.mat","df")
    df = df2;
    save(path2_+"data_2f50r.mat","df")
    df = df4;
    save(path2_+"data_2f50p.mat","df")
elseif number == 4
    df = spon;
    save(path2_+"data_2rsp.mat","df")
    df = df1;
    save(path2_+"data_2r0.mat","df")
    df = df3;
    save(path2_+"data_2r100p.mat","df")
    df = df5;
    save(path2_+"data_2r100r.mat","df")
    df = df2;
    save(path2_+"data_2r50r.mat","df")
    df = df4;
    save(path2_+"data_2r50p.mat","df")
    %df = df6;
    %save(path2_+"data_2r25r.mat","df")
elseif number == 5
    df = spon;
    save(path2_+"data_clsp.mat","df")
    df = df1;
    save(path2_+"data_cl0.mat","df")
    df = df3;
    save(path2_+"data_cl100p.mat","df")
    df = df5;
    save(path2_+"data_cl100r.mat","df")
    df = df2;
    save(path2_+"data_cl50r.mat","df")
    df = df4;
    save(path2_+"data_cl50p.mat","df")
    %df = df6;
    %save(path2_+"data_cl25r.mat","df")
elseif number == 6
    df1 = df;
    %df = rmfield(df,{'csp','cre'});
    df = rmfield(df,{'c50','c75','c100','c150','c200','c300','c400'});
    %save(path2_+"data_ex2.mat","df","-v7.3")
    save(path2_+"data_ex2_spon.mat","df","-v7.3")
    
    df = rmfield(df1,{'csp'});
    save(path2_+"data_ex2.mat","df","-v7.3")
end

