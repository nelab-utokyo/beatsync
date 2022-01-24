%% load music data
cnames = ["s75","s100","s200","s400","off"];

load('/Users/yoshiki/research/Tlab/code/sound/trigger/click_trg2.mat')
allbeats = reshape(onbeat.',1,33*4);
allbeats(isnan(allbeats)) = [];
allbeats2 = reshape(onbeat.',1,33*4);

spath = "/Users/yoshiki/research/Tlab/behavior/code/accelerator/software/sound/speed/";

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
%{
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
%}
%%

load("dno.mat")
allfiles = [dno1,dno2,dno3,dno4,dno5,dno6,dno7,dno8,dno9,dno11];

%% (jerk beat contrast)
mkinds = ["s75", "s100", "s200", "s400"];

indvdata = struct();

gap = -30:1:30;
onratio = 0.4;

for k = 1:size(allfiles,2)
    files = allfiles(1:3,k);
    prep = zeros(length(gap),4,2,length(files)*2);
    monamp = struct("s75",prep, "s100",prep, "s200",prep, "s400", prep);
    monamp2 = struct("s75",prep, "s100",prep, "s200",prep, "s400", prep);

    for j = 1:length(files)
        load(files(j).folder + "/" + files(j).name)
        dt = mean(diff(stamp));
        disp(files(j).name)
        disp([k,j,1/dt])
        dt2 = 1/40;
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
        for i = 4:4  %acceleration dimensions
            for n = 1:size(monoff,1)  %trial num
                for q = 1:1 %1: correct beat, 2: randomized beat(not used)
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
                            tmp2_ = diff(acc(:,st2-1:en),1,2)/ dt;
                            %tmp2_ = (zscore(acc(:,st2:monoff(n,3)),0,2));
                            tmp2 = sqrt(tmp2_(1,:).^2+tmp2_(2,:).^2+tmp2_(3,:).^2);
                        end                        
                        onb = mean(tmp2(tmp==1));
                        offb = mean(tmp2(tmp==0));
                        if q == 1
                            monamp.(mkinds(monoff(n,4)))(g,i,:,2*j-1+1*(n>4)) = reshape([onb,offb],[1,1,2]); 
                        elseif q > 1
                            monamp2.(mkinds(monoff(n,4)))(g,i,:,2*j-1+1*(n>4)) =  reshape([onb,offb],[1,1,2]);
                        end
                    end
                end
            end
        end
    end
    indvdata.("no"+k) = monamp;
    indvdata.("no2"+k) = monamp2;
end


%% day difference

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.2, 1]);
lnames = ["75%","100%","200%","400%"];

row = 4; % subplot行数
col = 1; % subplot列数
left_m = 0.15; % 左側余白の割合
bot_m = 0.05; % 下側余白の割合
ver_r = 1.2; % 縦方向余裕 (値が大きいほど各axes間の余白が大きくなる)
col_r = 1.15; % 横方向余裕 (値が大きいほど各axes間の余白が大きくなる)

for n = 1:4
   cond = n;
   ax(n) = axes('Position',...
      [(1-left_m)*(mod(n-1,col))/col + left_m ,...
      (1-bot_m)*(1-ceil(n/col)/(row)) + bot_m ,...
      (1-left_m)/(col*col_r ),...
      (1-bot_m)/(row*ver_r)]...
      );
    odr = [];
    allind = zeros(10,3);
    for k2 = 1:3
        for k = 1:10
            %tmp_ = squeeze(mean(indvdata.("no"+k).(mkinds(cond))(3:5,4,k2*2-1:k2*2,:)));
            tmp = squeeze(mean(indvdata.("no"+k).(mkinds(cond))(31,4,:,k2*2-1:k2*2),4));
            allind(k,k2) = (tmp(1) - tmp(2))./(tmp(1) + tmp(2));
        end
    end
    disp(["kruskalwallis",kruskalwallis(allind, [], "off")]) %ns ns ns ns
    %disp(friedman(allind, 1, "off"))       %ns ns ns ns 
    %p = signrank(allind(:,1),allind(:,3));
    p = ranksum(allind(:,1),allind(:,2),"tail","right");
    bx = boxplot(allind, "Widths",0.7);
    set(bx,'LineWidth',2)
    hold on
    for k = 1:10
        scatter((1:3)+rand(1,3)*0.2-0.1 ,allind(k,:),25,"og", "filled")
        hold on
    end
    disp([cond,mean(allind,1),p])  %0.089, 0.0452, 0.0757, 0.3847
    X = [ones(1,3); 1:3]';
    Y = X\mean(allind)';
    %plot(1:132, mean(allind),"g","Linewidth",1.5)
    %plot(1:3, X*Y, "r","Linewidth",3)
    plot([0,4], [0,0], "k--","Linewidth",1)
    
    cnt = 0;
    cand = [1,2];
    xt = [1,2,3];
    yt = 0.17;
    i = 1;
    plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt*(0.98+cnt/20), '-k', 'LineWidth',1.7)
    sign = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
    if sign >= 1
        for s = 1:sign
            plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2 - 0.5)/4, yt*(0.98+cnt/20) + 0.03, '*k',"markersize",18, 'LineWidth',1.7)
        end
    else
        text(mean(xt([cand(i,1) cand(i,2)]))-0.25, yt*(0.98+cnt/20)+0.045, "n.s.","Fontsize",25)
    end
    plot([1;1]*cand(i,1),[0.9+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
    plot([1;1]*cand(i,2),[0.9+cnt/20,0.98+cnt/20]*yt,'-k', 'LineWidth',1.7)
    ylim([-0.1,0.25])
    xlim([0.5,3.5])
    yticks(-0.2:0.1:0.2)
    xticks(1:3)
    xticklabels(["Day 1","Day 2","Day 3"])
    ax = gca;
    ax.YAxis.FontSize = 15;
    ax.XAxis.FontSize = 22;
    title(lnames(cond),"Fontsize",25)
    %xlabel("Beat index","Fontsize",20)
    %if cond == 1
    %    ylabel("Mean jerk difference","Fontsize",23)
    %end
end

%saveas(gca,"./summary3/beat_day.fig")
%saveas(gca,"./summary3/beat_day.png")

%% all individual data checker

names = ["x", "y", "z", "a"];
lnames = ["75%","100%","200%","400%"];

count = 0;
h = [];
knd = 2;
dim = 4;

chcker = zeros(length(gap),size(allfiles,2),length(files)*2);

for k = 1:size(allfiles,2)
    for i = 1:length(files)*2
        tmp = squeeze(indvdata.("no"+k).(mkinds(knd))(:,dim,:,i));
        chcker(:,k,i) = (tmp(:,1) - tmp(:,2))./(tmp(:,1) + tmp(:,2));
    end
end
figure()
%imagesc(squeeze(max(chcker)))
imagesc(squeeze(chcker(31,:,:)))


%% synchrony shift per trial

figure()
for knd = 1:4
    chcker = zeros(length(gap),size(allfiles,2),length(files)*2);
    for k = 1:size(allfiles,2)
        for i = 1:length(files)*2
            tmp = squeeze(indvdata.("no"+k).(mkinds(knd))(:,dim,:,i));
            chcker(:,k,i) = (tmp(:,1) - tmp(:,2))./(tmp(:,1) + tmp(:,2));
        end
    end
    tmp = squeeze(max(chcker));
    %tmp = squeeze((chcker(31,:,:)));
    errorbar(1:length(files)*2,mean(tmp,1),std(tmp,[],1)/sqrt(size(allfiles,2)))
    hold on
end


%% whole

ani = 6;
trial = 2;

load(allfiles(ceil(trial/2), ani).folder + "/" + allfiles(ceil(trial/2), ani).name)

exs = find(ontime(:,3) == 1);
exs = exs(trial);

t1 = ontime(exs,1) - ontime(exs,2)/1000;
t2 = offtime(exs,1) - (mlns1(offtime(exs,3)+1) - offtime(exs,2))/1000;


[M,I1]=min(abs(stamp - t1));
[M,I2]=min(abs(stamp - t2));

time = stamp(I1:I2-1)- stamp(I1);


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.5]);


clck = zeros(size(time))-0.2;
for i = 1:132
    [M,I] = min(abs(time  - aclick.s100(i)/1000));
    clck(I) = 10;
    plot([1,1]*time(I), [-1,10], "k")
    hold on
end

hs = [];
for i = 1:3
    h = plot(time,abs(diff(acc(i,I1:I2)))+2 - i/2, "linewidth",1.5);
    hs = [hs,h];
end

%legend(hs, ["|da_x|","|da_y|","|da_z|"], "Fontsize",20)
legend(hs, ["X","Y","Z"], "Fontsize",25, "location", "northeastoutside")

%set(gca, "Xticklabels",[],"Yticklabels",[])
ylim([0.3,2.5])

%% zoom in rat6 (beat synchrony)

clrs = [1,0.5,0;0.5,1,0;0,0.5,1]; % with color

%clrs = [0,0,0;0.35,0.35,0.35;0.7,0.7,0.7];  % black or white

ani = 6;
trial = 2;

load(allfiles(ceil(trial/2), ani).folder + "/" + allfiles(ceil(trial/2), ani).name)

exs = find(ontime(:,3) == 1);
exs = exs(trial);

t1 = ontime(exs,1) - ontime(exs,2)/1000;
t2 = offtime(exs,1) - (mlns1(offtime(exs,3)+1) - offtime(exs,2))/1000;

[M,I1]=min(abs(stamp - t1));
[M,I2]=min(abs(stamp - t2));

time = stamp(I1:I2-1)- stamp(I1);


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.5]);


clck = zeros(size(time))-0.2;
for i = 1:132
    [M,I] = min(abs(time  - aclick.s100(i)/1000));
    clck(I) = 10;
    plot([1,1]*time(I), [-1,10], "k")
    hold on
end

I3 = 34037;
I4 = 34455;

hs = [];
accgp = [0.4,0.75,2];
accgp2 = [-0.8,0.55,0.8];
for i = 1:3
    %h = plot(time((I3:(I4-1))-I1),abs(diff(acc(i,I3:I4)))+2 - i/2, "linewidth",1.5, "color", clrs(i,:));
    %h2.Color = [h2.Color 0.5];
    %h = plot(time((I3-1:(I4-1))-I1),abs((acc(i,I3:I4)))+2 - accgp(i), "linewidth",1.5, "color", clrs(i,:));
    h = plot(time((I3-1:(I4-1))-I1),acc(i,I3:I4)+2 - accgp2(i), "linewidth",1.5, "color", clrs(i,:));
    hs = [hs,h];
end

%legend(hs, ["|da_x|","|da_y|","|da_z|"], "Fontsize",20)
legend(hs, ["X","Y","Z"], "Fontsize",25, "location", "northeastoutside")

ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ylim([0.3,2.5])
ylim([0.2,3.3])
%ylim([0,4])
xlim([time(I3-I1),time(I4-I1)])
%whitebg("white")

%% zoom in for video rat 6 (beat synchrony)

figpath = "/Users/yoshiki/research/Tlab/behavior/acceleration/20201207/movies/figures/";

clrs = [1,0.5,0;0.5,1,0;0,0.5,1];
ani = 6;
trial = 2;

load(allfiles(ceil(trial/2), ani).folder + "/" + allfiles(ceil(trial/2), ani).name)

exs = find(ontime(:,3) == 1);
exs = exs(trial);

t1 = ontime(exs,1) - ontime(exs,2)/1000;
t2 = offtime(exs,1) - (mlns1(offtime(exs,3)+1) - offtime(exs,2))/1000;

[M,I1]=min(abs(stamp - t1));
[M,I2]=min(abs(stamp - t2));

time = stamp(I1:I2-1)- stamp(I1);
I3_ = 34037;
I4_ = 34395;

%gp = (I4_ - I3_)/200;
gp = 2;
gtime = (0:length(time)-1);

for j = 1:1:1 %179
    now = gp*j;
    %figure("Visible","off")
    fig = figure("Visible","on");
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.25, 0.25]);

    I5 = round(I3_ + now);
    I3 = I5 - 150;
    I4 = I5 + 150;
    hs = [];
    for i = 1:3
        h = plot(-149:150,abs(diff(acc(i,I3:I4)))+2 - i/2, "linewidth",2, "color", clrs(i,:));
        hold on
        hs = [hs,h];
    end
    plot([1,1]*0,[0,5], "r", "linewidth",4)
    %legend(hs, ["|da_x|","|da_y|","|da_z|"], "Fontsize",20)
    legend(hs, ["X","Y","Z"], "Fontsize",25, "location", "northeastoutside", "box","off")
    set(gcf,'InvertHardcopy','off')
    set(gca,'Color','none');

    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ylim([0.3,2.5])
    xlim([-150,150])
    rectangle("pos",[-160,1,1,10],"Facecolor","k")
    whitebg(fig,[0.3 .3 .3])
    %print(figpath + "plot_"+j, "-dpng")
    %close
end



%% zoom in rat4

clrs = [1,0.5,0;0.5,1,0;0,0.5,1];
ani = 6;
trial = 1;

load(allfiles(ceil(trial/2), ani).folder + "/" + allfiles(ceil(trial/2), ani).name)

exs = find(ontime(:,3) == 1);
exs = exs(trial);

t1 = ontime(exs,1) - ontime(exs,2)/1000;
t2 = offtime(exs,1) - (mlns1(offtime(exs,3)+1) - offtime(exs,2))/1000;

[M,I1]=min(abs(stamp - t1));
[M,I2]=min(abs(stamp - t2));

time = stamp(I1:I2-1)- stamp(I1);
gtime = (0:length(time)-1);


figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.6]);

I3 = I1 + 1800;
%I3 = I1+2200;
I4 = I1 + 5100;


hs = [];
dif = [0,2,3];
for i = 1:3
    h = plot((I3:(I4))-I1,((acc(i,I3:I4)))+2 - dif(i), "linewidth",1.5, "color", clrs(i,:));
    hold on
    hs = [hs,h];
end
plot([1,1]*0, [-1,10], "k","linewidth",3)
plot([1,1]*(I2-I1), [-10,10], "k","linewidth",5)

%legend(hs, ["|da_x|","|da_y|","|da_z|"], "Fontsize",20)
legend(hs,["X","Y","Z"], "Fontsize",25, "location", "northeastoutside")

ax = gca;

%ax.YAxis.Visible = 'off';
ylim([-2.5,3.4])
xlim([I3-I1,I4-I1])


ff = 5*(I4-I1)/(stamp(I4) - stamp(I3));
xticks((I2-I1)-ff*5:ff:(I2-I1)+ff*6)
xticklabels(-25:5:30)
%ax.XAxis.Visible = 'off';


%%

dt = (stamp(I4) - stamp(I3))/(I4-I3);
%dt = 1/43;
range = [0.41,19];
[wt,frq] = cwt(acc(1,:), 1/dt,"FrequencyLimits",range);

frqln = 7;
gap = 8;
tmp = [];
tmp_ = abs(wt(:,I3:I4));

for l = 1:frqln
    tmp = [tmp;mean(tmp_((gap*l-gap+1):gap*l,:))];
end

figure()

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.6]);

imagesc(tmp_)
colormap default
hold on
plot([1,1]*(I2-I3), [-1,length(frq)+1], "k","linewidth",3)

yticks(4:gap:55)
yticklabels(round(frq(4:gap:55)*10)/10)
ff = round(10/dt);
xticks((I2-I3)-ff*5:ff:(I2-I3)+ff*6)
xticklabels(-50:10:60)

