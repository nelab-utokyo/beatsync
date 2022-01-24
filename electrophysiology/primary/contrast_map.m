%% 
%%%%%%%%%%%%%%%%%%
% rhythmic click
%%%%%%%%%%%%%%%%%%

paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006"];

labels = [60,120,240,480] + " BPM";

ctrf_ = [];
aaf = [];
a1 = [];
cb = [];
cfs = [];
lats = [];
neunum = [0];

for p = 1:length(paths)
    load("./"+paths(p)+"/buildmat/responsep.mat")
    load("./"+paths(p)+"/buildmat/responsef.mat")
    load("./"+paths(p)+"/shuf_param.mat")
    load("./"+paths(p)+"/CF_latency.mat")
    ffs(rmv,:,:) = [];
    TC(:,rmv) = [];
    ctrf = [];
    for n = 1:size(ffs,1)
        ctrf = [ctrf;(ffs(n,:,1)-ffs(n,:,2))./(ffs(n,:,1)+ffs(n,:,2))];
    end
    ctrf_ = [ctrf_;ctrf];
    cb = [cb; corebelt(TC)];
    aaf = [aaf; (a1aaf(TC)==-1)];
    a1 = [a1; (a1aaf(TC)==1)];
    cfs = [cfs; CF(TC)];
    lats = [lats; latency(TC)];
    neunum = [neunum, size(ctrf_,1)];
end


%%

temp1 = [];
for p = 1:length(paths)
    pkup = ctrf_(neunum(p)+1:neunum(p+1),:);
    temp1 = [temp1;nanmean(pkup(find(cb(neunum(p)+1:neunum(p+1))==1),:))];
end
%temp1 = ctrf_(find(aaf==1),:);

temp2 = [];
for p = 1:length(paths)
    pkup = ctrf_(neunum(p)+1:neunum(p+1),:);
    temp2 = [temp2;nanmean(pkup(find(cb(neunum(p)+1:neunum(p+1))==0),:))];
end
%temp2 = ctrf_(find(a1==1),:);


cand = [1,2;];
ps = ones(length(labels),length(cand));
for i = 1:length(labels)
    p = anova1([temp1(:,i),temp2(:,i)],[],"off");
    if p < 0.05
        [p1,h] = ranksum(temp1(:,i),temp2(:,i));
        %[p1,h] = signrank(temp1(:,i),temp2(:,i));
        ps(i,1) = p1; 
        disp([i,p,p1])
    else
        disp([i,p])
    end
end

%% 2 way anova for region and BPM

anv2 = [temp1;temp2];
disp(anova2(anv2,length(paths),"off"))  %2.3989e-10, 0.1756
[~,~,stats] = anova2(anv2,length(paths),"off");
c = multcompare(stats);
disp(c(:,[1,2,end]))

%%
cand_ = [1,2;2,3;3,4;];
disp(anova1(temp1,[],"off"))  %8.1074e-07
disp(anova1(temp2,[],"off"))  %1.3356e-04

ps2 = zeros(2,length(cand_));
for i = 1:length(cand_)
    ps2(1,i) = signrank(temp1(:,cand_(i,1)),temp1(:,cand_(i,2)));   %0.0039, 0.4258, 0.0039
    ps2(2,i) = signrank(temp2(:,cand_(i,1)),temp2(:,cand_(i,2)));   %0.0547, 0.0273, 0.0117
    %ps2(1,i) = ranksum(temp1(:,cand_(i,1)),temp1(:,cand_(i,2)));
    %ps2(2,i) = ranksum(temp2(:,cand_(i,1)),temp2(:,cand_(i,2)));
    disp([cand_(i,:),ps2(:,i)'])
end

%%
ttls = ["Core","Belt"];
labels = [60,120,240,480];
clrs = [1,0.3,0.3;0.3,0.3,1];
%clrs = [1,0.6,0.6;0.1,0.1,0.5];
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);
pts = [];
for i = 1:2
    if i == 1
        temp = temp1;
    elseif i == 2
        temp = temp2;
    else
        temp = temp3;
    end
    h = boxplot(temp,"Colors",clrs(i,:),"Widths",0.2,"Positions",(1:4)+(i/4-0.375),"OutlierSize",3,"Symbol","k+");
    set(h,'LineWidth',2)
    hold on
    for j = 1:size(temp,1)
        scatter((1:4)+(i/4-0.375) + (rand(1,4)-0.5)*0.05, temp(j,:), 50, "k", "filled")
    end
    b = plot([-1,-1],[0,1],"Color",clrs(i,:),"linewidth",1.7);
    pts = [pts,b];
end


yt = [0.44,0.39];
xt_ = 1:6;
for l = size(ps2,1):-1:1
    xt = xt_ + (l - 1.5)/4;
    cnt = 0;
    cnt2 = 0;
    for i = 1:length(cand_)
        if ps2(l,i) < 0.05
            plot(xt([cand_(i,1) cand_(i,2)]), [1 1]*yt(l)+cnt, '-',"color",clrs(l,:), 'LineWidth',2)
            sign = (1 + 1*(ps2(l,i) < 0.01) + 1*(ps2(l,i) < 0.001));
            for s = 1:sign
                plot(mean(xt([cand_(i,1) cand_(i,2)]))+(s-sign/2-0.5)/7, yt(l)+cnt+0.01, '*k',"markersize",15, 'LineWidth',1.7)
            end
            plot([1;1]*xt(cand_(i,1)),[-0.02,0]+yt(l)+cnt,'-',"color",clrs(l,:), 'LineWidth',2)
            plot([1;1]*xt(cand_(i,2)),[-0.02,0]+yt(l)+cnt,'-',"color",clrs(l,:), 'LineWidth',2)
            cnt = cnt + cnt2;
        else
            dy = 0.0;
            plot(xt([cand_(i,1) cand_(i,2)]), [1 1]*yt(l)+cnt-dy, '-',"color",clrs(l,:), 'LineWidth',2)
            sign = (1 + 1*(ps2(l,i) < 0.01) + 1*(ps2(l,i) < 0.001));
            text(mean(xt([cand_(i,1) cand_(i,2)]))-0.1, yt(l)+cnt+0.018-dy, "n.s.","Fontsize",25)
            plot([1;1]*xt(cand_(i,1)),[-0.02,0]+yt(l)+cnt-dy,'-',"color",clrs(l,:), 'LineWidth',2)
            plot([1;1]*xt(cand_(i,2)),[-0.02,0]+yt(l)+cnt-dy,'-',"color",clrs(l,:), 'LineWidth',2)
        end
    end
end

%{
yt = [0.26,0.38,0.33,0.2];
xt_ = [-0.25,0,0.25];
for l = 1:size(ps,1)
    xt = xt_ + l;
    cnt = 0;
    cnt2 = 0.018;
    for i = 1:length(cand)
        if ps(l,i) < 0.05
            plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt(l)+cnt, '-k', 'LineWidth',1.5)
            sign = (1 + 1*(ps(l,i) < 0.01) + 1*(ps(l,i) < 0.001));
            for s = 1:sign
                plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2-0.5)/10, yt(l)+cnt+0.008, '*k',"markersize",10, 'LineWidth',1.5)
            end
            plot([1;1]*xt(cand(i,1)),[-0.01,0]+yt(l)+cnt,'-k', 'LineWidth',1.5)
            plot([1;1]*xt(cand(i,2)),[-0.01,0]+yt(l)+cnt,'-k', 'LineWidth',1.5)
            cnt = cnt + cnt2;
        end
    end
end
%}

xticks(1:4)
xticklabels(labels)
yticks(-0.2:0.1:0.5)

ylim([-0.12,0.48])
xlim([0.3,4.7])
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;
%xlabel("BPM", "Fontsize",28)
%ylabel("Beat contrast", "Fontsize",28)
legend(pts,ttls,"Fontsize",28,"Location","west")

%saveas(gcf,"./summary/beat_contrast_f_region2.fig")
%saveas(gcf,"./summary/beat_contrast_f_region2.png")

%saveas(gcf,"./summary/beat_contrast_f_region2_bw.fig")
%saveas(gcf,"./summary/beat_contrast_f_region2_bw.png")
%%
%%%%%%%%%%%%%%%%%%
% music (calculate musical on beat response)
%%%%%%%%%%%%%%%%%%


paths = ["20200908","20200910","20200915","20200917","20201006","20210106","20210107"];
%paths = ["20200908","20200910","20200915","20200917","20201006","20201220","20201222","20210106","20210107","20210108_1","20210108_2"];
%paths = ["20200908","20200910","20200915","20200917","20201006","20201220","20210106","20210107","20210108_1"];

cands = [75,100,200,300,400];

%paths = ["20201220","20201222","20210106","20210107","20210108_1","20210108_2"];
%cands = [50,75,100,150,200,400];


labels = cands + " BPM";

names = cands +"%(" + round((cands/100)*132) + ")";
bresall = [];
nbresall = [];

aaf = [];
a1 = [];
cb = [];
cfs = [];
lats = [];
neunum = [0];

for p = 1:length(paths)
    load("./"+paths(p)+"/shuf_param.mat")
    load("./"+paths(p)+"/CF_latency.mat")
    load("./"+paths(p)+"/buildmat/responsem.mat")
    TC(rmv) = [];
    bres(:,rmv,:) = [];
    nbres(:,rmv,:) = [];
    bresall = [bresall; permute(bres,[2,1,3])];
    nbresall = [nbresall; permute(nbres,[2,1,3])];
    neunum = [neunum, size(bresall,1)];
    cb = [cb; corebelt(TC)];
    cfs = [cfs; CF(TC)];
    aaf = [aaf; (a1aaf(TC)==-1)];
    a1 = [a1; (a1aaf(TC)==1)];
    lats = [lats; latency(TC)];
    disp(length(cb))
end

contrast = (bresall(:,:,5) - nbresall(:,:,5))./(bresall(:,:,5) + nbresall(:,:,5));


%%

temp1 = [];
for p = 1:length(paths)
    pkup = contrast(neunum(p)+1:neunum(p+1),:);
    temp1 = [temp1;nanmean(pkup(find(cb(neunum(p)+1:neunum(p+1))==1),:))];
    %temp1 = [temp1;nanmean(pkup(find(cfs(neunum(p)+1:neunum(p+1))>12),:))];
end
%temp1 = contrast(find(aaf==1),:);

temp2 = [];
for p = 1:length(paths)
    pkup = contrast(neunum(p)+1:neunum(p+1),:);
    temp2 = [temp2;nanmean(pkup(find(cb(neunum(p)+1:neunum(p+1))==0),:))];
    %temp2 = [temp2;nanmean(pkup(find(cfs(neunum(p)+1:neunum(p+1))<6),:))];
end
%temp2 = contrast(find(a1==1),:);


cand = [1,2;];
ps = zeros(length(labels),length(cand));
for i = 1:length(labels)
    p = anova1([temp1(:,i),temp2(:,i)],[],"off");
    if p < 0.05
        [p1,h] = ranksum(temp1(:,i),temp2(:,i));
        ps(i,:) = p1; 
        disp([i,p,p1])
    end
end


%% 2way anova for region and BPM

anv2 = [temp1;temp2];
%anv2 = [temp3;temp4;];
disp(anova2(anv2,length(paths),"off"))   %7.712e-5, 5.9974e-18
anova2(anv2,length(paths),"off")
[~,~,stats] = anova2(anv2,length(paths),"off");
c1 = multcompare(stats,"Estimate","column");
c2 = multcompare(stats,"Estimate","row");

disp(c1(:,[1,2,end]))
disp(c2(:,[1,2,end]))

%%
%cand_ = [4,5;3,4;2,3;1,2;2,5;];
cand_ = [1,2;2,3;3,4;4,5];
disp(anova1(temp1,[],"off"))  %5.9974e-18
disp(anova1(temp2,[],"off"))  %2.4551e-06
ps2 = zeros(2,length(cand_));

for i = 1:length(cand_)    
    ps2(1,i) = signrank(temp1(:,cand_(i,1)),temp1(:,cand_(i,2)));    %0.2188, 0.0156, 0.0156, 0.0156
    ps2(2,i) = signrank(temp2(:,cand_(i,1)),temp2(:,cand_(i,2)));    %0.1562, 0.0156, 0.0312, 0.1562
    disp([cand_(i,:),ps2(:,i)'])
end

%%
ttls = ["Core","Belt"];
clrs = [1,0.3,0.3;0.3,0.3,1];
clrs = [1,0.6,0.6;0.1,0.1,0.5];

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);
pts = [];
for i = 1:2
    if i == 1
        temp = temp1;
    elseif i == 2
        temp = temp2;
    else
        temp = temp3;
    end
    h = boxplot(temp,"Colors",clrs(i,:),"Widths",0.2,"Positions",(1:length(labels))+(i/4-0.375),"OutlierSize",3,"Symbol","k+");
    set(h,'LineWidth',1.7)
    hold on
    for j = 1:size(temp,1)
        scatter((1:5)+(i/4-0.375) + (rand(1,5)-0.5)*0.05, temp(j,:), 50, "k", "filled")
    end
    b = plot([-1,-1],[0,1],"Color",clrs(i,:),"linewidth",1.7);
    pts = [pts,b];
end


yt = [0.52,0.48];
xt_ = 1:length(labels);
for l = size(ps2,1):-1:1
    xt = xt_ + (l - 1.5)/4;
    cnt = 0;
    cnt2 = 0.025;
    for i = 1:length(cand_)
        if ps2(l,i) < 0.05
            cnt = cnt + cnt2*(i == 5);
            plot(xt([cand_(i,1) cand_(i,2)]), [1 1]*yt(l)+cnt, '-',"color",clrs(l,:), 'LineWidth',2)
            sign = (1 + 1*(ps2(l,i) < 0.01) + 1*(ps2(l,i) < 0.001));
            for s = 1:sign
                plot(mean(xt([cand_(i,1) cand_(i,2)]))+(s-sign/2-0.5)/10, yt(l)+cnt+0.01, '*k',"markersize",15, 'LineWidth',1.7)
            end
            plot([1;1]*xt(cand_(i,1)),[-0.01,0]+yt(l)+cnt,'-',"color",clrs(l,:), 'LineWidth',2)
            plot([1;1]*xt(cand_(i,2)),[-0.01,0]+yt(l)+cnt,'-',"color",clrs(l,:), 'LineWidth',2)
        else
            dy = 0.0;
            plot(xt([cand_(i,1) cand_(i,2)]), [1 1]*yt(l)+cnt-dy, '-',"color",clrs(l,:), 'LineWidth',2)
            sign = (1 + 1*(ps2(l,i) < 0.01) + 1*(ps2(l,i) < 0.001));
            text(mean(xt([cand_(i,1) cand_(i,2)]))-0.1, yt(l)+cnt+0.018-dy, "n.s.","Fontsize",25)
            plot([1;1]*xt(cand_(i,1)),[-0.01,0]+yt(l)+cnt-dy,'-',"color",clrs(l,:), 'LineWidth',2)
            plot([1;1]*xt(cand_(i,2)),[-0.01,0]+yt(l)+cnt-dy,'-',"color",clrs(l,:), 'LineWidth',2)
        end
    end
end


%{
yt = [0.5,0.55,0.61,0.61,0.5];
xt_ = [-0.25,0,0.25];
for l = 1:size(ps,1)
    xt = xt_ + l;
    cnt = 0;
    cnt2 = 0.025;
    for i = 1:length(cand)
        if ps(l,i) < 0.05
            plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt(l)+cnt, '-k', 'LineWidth',1.5)
            sign = (1 + 1*(ps(l,i) < 0.01) + 1*(ps(l,i) < 0.001));
            for s = 1:sign
                plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2-0.5)/10, yt(l)+cnt+0.008, '*k',"markersize",10, 'LineWidth',1.3)
            end
            plot([1;1]*xt(cand(i,1)),[-0.01,0]+yt(l)+cnt,'-k', 'LineWidth',1.5)
            plot([1;1]*xt(cand(i,2)),[-0.01,0]+yt(l)+cnt,'-k', 'LineWidth',1.5)
            cnt = cnt + cnt2;
        end
    end
end
%}

xticks(1:length(cands))
xticklabels(cands + "%")
yticks(0.1:0.1:0.8)

ylim([0.08,0.56])
xlim([0.3,length(labels) + 0.7])
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;
%xlabel("Play speed(BPM)", "Fontsize",32)
%ylabel("Beat contrast", "Fontsize",32)
legend(pts,ttls,"Fontsize",28,"Location","Southwest")

%saveas(gcf,"./summary/beat_contrast_region2.fig")
%saveas(gcf,"./summary/beat_contrast_region2.png")

%saveas(gcf,"./summary/beat_contrast_region2_bw.fig")
%saveas(gcf,"./summary/beat_contrast_region2_bw.png")


%%
ttls = ["AAF","A1","Belt"];
clrs = [1,0.3,0.3;0.3,0.6,0.3;0.3,0.3,1];
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.6]);
pts = [];
for i = 1:3
    if i == 1
        temp = temp1;
    elseif i == 2
        temp = temp2;
    else
        temp = temp3;
    end
    boxplot(temp,"Colors",clrs(i,:),"Widths",0.2,"Positions",(1:5)+(i/4-0.5),"OutlierSize",3,"Symbol","k+");
    hold on
    b = plot([-1,-1],[0,1],"Color",clrs(i,:),"linewidth",1.5);
    pts = [pts,b];
end

yt = [0.5,0.55,0.55,0.52,0.4];
xt_ = [-0.25,0,0.25];
for l = 1:size(ps,1)
    xt = xt_ + l;
    cnt = 0;
    cnt2 = 0.025;
    for i = 1:length(cand)
        if ps(l,i) < 0.05
            plot(xt([cand(i,1) cand(i,2)]), [1 1]*yt(l)+cnt, '-k', 'LineWidth',1.5)
            sign = (1 + 1*(ps(l,i) < 0.01) + 1*(ps(l,i) < 0.001));
            for s = 1:sign
                plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2-0.5)/10, yt(l)+cnt+0.008, '*k',"markersize",10, 'LineWidth',1.3)
            end
            plot([1;1]*xt(cand(i,1)),[-0.01,0]+yt(l)+cnt,'-k', 'LineWidth',1.5)
            plot([1;1]*xt(cand(i,2)),[-0.01,0]+yt(l)+cnt,'-k', 'LineWidth',1.5)
            cnt = cnt + cnt2;
        end
    end
end
%}

xticks(1:6)
xticklabels(names)
yticks(0.1:0.1:0.8)

ylim([0.13,0.62])
xlim([0.3,5.7])
ax = gca;
ax.YAxis.FontSize = 21;
ax.XAxis.FontSize = 21;
xlabel("Play speed(BPM)", "Fontsize",25)
ylabel("Beat contrast", "Fontsize",25)
legend(pts,ttls,"Fontsize",20,"Location","Northeast")

%saveas(gcf,"./summary/beat_contrast_region2x.fig")
%saveas(gcf,"./summary/beat_contrast_region2x.png")

