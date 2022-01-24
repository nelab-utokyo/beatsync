%%
%%%%%%%%%%%%
% calculatebeat contrast of p & f 
%%%%%%%%%%%%

paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006"];

labels = [60,120,240,480];
splfp = 1;
ctrf = [];
resf = [];
resp = [];
neunum = [0];
spns = [];
fails = [];
rmvs = [];
for p = 1:length(paths)
    load("./"+paths(p)+"/buildmat/responsep.mat")
    load("./"+paths(p)+"/buildmat/responsef.mat")
    ffs(rmv,:,:) = [];
    diffs(rmv,:,:) = [];
    failrs(rmv,:,:) = [];
    resf = [resf; ffs];
    resp = [resp; diffs];
    fails = [fails; failrs];
    rmvs = [rmvs;rmv];
    for n = 1:size(ffs,1)
        ctrf = [ctrf;(ffs(n,:,1)-ffs(n,:,2))./(ffs(n,:,1)+ffs(n,:,2))];
    end
    disp(size(ctrf))
    neunum = [neunum, size(ctrf,1)];
end

%% calculate the significance

temp = [];
for i = 1:length(neunum)-1
    temp = [temp; nanmean(ctrf(neunum(i)+1:neunum(i+1),:))];
end
cand = [1,2;1,3;1,4;2,3;2,4;3,4];
sign = zeros(length(cand),1);
ps = zeros(length(cand),1);
for i = 1:length(cand)
    [p,h] = signrank(temp(:,cand(i,1)),temp(:,cand(i,2)));
    ps(i) = p;
    disp([cand(i,:),p])
    sign(i) = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
end

%% plot beat contrast of f
%%%%%%%%%%%%%%%%
% calculate rhythmic click beat contrast
%%%%%%%%%%%%%%%%
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);
temp = [];
for i = 1:length(neunum)-1
    temp = [temp; nanmean(ctrf(neunum(i)+1:neunum(i+1),:))];
end

disp("Anova = " + anova1(temp,labels,"off"))  %4.0816e-07
%errorbar(1:length(labels),nanmean(ctrf),nanstd(ctrf)/sqrt(size(ctrf,1)),"o-", 'LineWidth',2)
h = boxplot(temp);
hold on
for i = 1:length(neunum)-1
    scatter((1:size(temp,2))+(rand(1,4)-0.5)*0,temp(i,:),100,"filled")
end
set(h,'LineWidth',2)
xlim([0.6,4.4])
xticks(1:length(labels))
xticklabels(labels)
yticks(-0.1:0.1:0.5)
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;

xt = get(gca, 'XTick');
yt = [0,0.33];
cnt = 0.5;
cnt2 = 2.5;
ylim([-0.1,0.35])


saveas(gca,"./summary/beat_contrast_f.fig")
saveas(gca,"./summary/beat_contrast_f.png")


%%
%%%%%%%%%%%%%%%%
% calculate musical on beat response
%%%%%%%%%%%%%%%%

paths = ["20200908","20200910","20200915","20200917","20201006","20210106","20210107"]; 
cands = [75,100,200,300,400];

names = cands +"%(" + round((cands/100)*132) + ")";

ctrf = [];
bresall = [];
nbresall = [];
neunum = [0];

for p = 1:length(paths)
    load("./"+paths(p)+"/shuf_param.mat")
    load("./"+paths(p)+"/buildmat/responsem.mat")
    disp([length(TC),size(bres,2)])
    bres(:,rmv,:) = [];
    nbres(:,rmv,:) = [];
    bresall = [bresall; permute(bres,[2,1,3])];
    nbresall = [nbresall; permute(nbres,[2,1,3])];
    neunum = [neunum, size(bresall,1)];
end

contrast = (bresall(:,:,5) - nbresall(:,:,5))./(bresall(:,:,5) + nbresall(:,:,5));


%% calculate the significance

temp = [];
for i = 1:length(paths)
    temp = [temp; nanmean(contrast(neunum(i)+1:neunum(i+1),:),1)];
end
disp("Anova = " + anova1(temp,names,"off"))  %3.6561e-05

ps = zeros(size(contrast,2))*NaN(1);
for i = 1:length(ps)
    for j = i+1:length(ps)
        [p,h] = signrank(temp(:, i),temp(:,j));
        ps(i,j) = p;
    end
end

disp("min = "+min(ps))
figure()

h = imagesc(ps);
hold on
colorbar()
set(gca,'ColorScale','log')
h.AlphaData = isfinite(ps); % NaNやInfを透明にする
h.Parent.Color = 'k';
for i = 1:length(ps)
    for j = i+1:length(ps)
        if ps(i,j) < 0.005
            plot([j-0.2,j,j+0.2], [i,i,i],'*k', "MarkerSize", 5)
        elseif ps(i,j) < 0.01
            plot([j-0.1,j+0.1], [i,i],'*k', "MarkerSize", 5)
        elseif ps(i,j) < 0.05
            plot(j, i,'*k', "MarkerSize", 5)
        end
    end
end
xticks(1:length(names))
yticks(1:length(names))
xticklabels(names)
yticklabels(names)
title("p value", "fontsize", 20)

%saveas(gca, "./summary/music_pval.png")
%saveas(gca, "./summary/music_pval.fig")

%% plot musical on-beat response of all position

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);

h = boxplot(temp);
set(h,'LineWidth',2)
hold on
for i = 1:size(temp,1)
    s1 = scatter(1:size(temp,2), temp(i,:), 100, "o","filled");
end

xt = get(gca, 'XTick');
yt = [0,0.42];
cnt = 0;
cnt2 = 0.0;
cand = [4,5;3,4;2,3;1,2];
%cand = [3,4;2,3;2,4;1,4;1,2;];
%
for i = 1:length(cand)
    p = ps(cand(i,1),cand(i,2));
    sign = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
    cnt = cnt + cnt2;
    if sign  > 0
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*max(yt)*0.98+cnt, '-k', 'LineWidth',2)
        for s = 1:sign
            plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2-0.4)/8, max(yt)*0.995+cnt, '*k',"markersize",15, 'LineWidth',1.7)
        end
        plot([1;1]*cand(i,1),[0.95,0.98]*max(yt)+cnt,'-k', 'LineWidth',2)
        plot([1;1]*cand(i,2),[0.95,0.98]*max(yt)+cnt,'-k', 'LineWidth',2)
    else
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*max(yt)*0.98+cnt, '-k', 'LineWidth',2)
        text(mean(xt([cand(i,1) cand(i,2)]))-0.1, max(yt)*0.98+0.01+cnt, "n.s.","Fontsize",25)
        plot([1;1]*cand(i,1),[0.95,0.98]*max(yt)+cnt,'-k', 'LineWidth',2)
        plot([1;1]*cand(i,2),[0.95,0.98]*max(yt)+cnt,'-k', 'LineWidth',2)
        cnt = cnt + cnt2;
    end
end

cand = [1,4;2,4];
yt = [0,0.21];
cnt = 0;
cnt2 = 0.018;

%
for i = 1:length(cand)
    p = ps(cand(i,1),cand(i,2));
    sign = 1*(p < 0.05) + 1*(p < 0.01) + 1*(p < 0.001);
    cnt = cnt + cnt2*(i >= 5);
    if sign  > 0
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*max(yt)+cnt, '-k', 'LineWidth',2)
        for s = 1:sign
            plot(mean(xt([cand(i,1) cand(i,2)]))+(s-sign/2-0.4)/8, max(yt)*0.97+cnt, '*k',"markersize",15, 'LineWidth',1.7)
        end
        plot([1;1]*cand(i,1),[1.0,1.05]*max(yt)+cnt,'-k', 'LineWidth',2)
        plot([1;1]*cand(i,2),[1.0,1.05]*max(yt)+cnt,'-k', 'LineWidth',2)
    else
        plot(xt([cand(i,1) cand(i,2)]), [1 1]*max(yt)+cnt, '-k', 'LineWidth',2)
        text(mean(xt([cand(i,1) cand(i,2)]))-0.1, max(yt)-0.008+cnt, "n.s.","Fontsize",25)
        plot([1;1]*cand(i,1),[1.0,1.05]*max(yt)+cnt,'-k', 'LineWidth',2)
        plot([1;1]*cand(i,2),[1.0,1.05]*max(yt)+cnt,'-k', 'LineWidth',2)
        cnt = cnt + cnt2;
    end
end

xticks(1:length(names))
xticklabels(cands + "%")
yticks(-0.4:0.1:1)
ylim([0.17,0.44])

xlim([0.5,length(names)+0.5])
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;

%ylabel("Beat contrast", "FontSize",28)
%xlabel("Play speed (BPM)", "FontSize",28)
%saveas(gca, "./summary/beat_contrast.png")
%saveas(gca, "./summary/beat_contrast.fig")


%% plot musical on-beat response of each position (with color)

names = cands + "%(" + round((cands/100)*132) + ")";
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);

data = zeros(length(neunum)-1,5,4);
for i = 1:length(neunum)-1
    st = neunum(i)+1;
    en = neunum(i+1);
    for j = 1:4
        temp = (bresall(st:en,:,j) - nbresall(st:en,:,5))./(bresall(st:en,:,j) + nbresall(st:en,:,5));
        data(i,:,j) = mean(temp);
    end
end

clors = [0.8500,0.3250,0.0980;0,0.4470,0.7410;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560];
mrk = ["o-","^-","d-","v-"];
h1s = [];
for i = 1:4
    temp = data(:,:,i);
    h1 = errorbar((1:length(cands))-0.25+0.1*i,nanmean(temp,1),nanstd(temp,0,1),mrk(i),"Color",clors(i,:),'LineWidth',2,"Markersize",20);
    set(h1, 'markerfacecolor', get(h1, 'color'))
    h1s = [h1s,h1];
    hold on
end

xticks(1:length(names))
xticklabels(cands + "%")
yticks(0:0.1:0.8)

xlim([0.5,length(names)+0.5])
ylim([0.08,0.5])
%legend(h1s,["pos.1","pos.2","pos.3","pos.4"], "location","south west","fontsize", 23);
[h,icons] = legend(h1s,["pos.1","pos.2","pos.3","pos.4"], "location","south west","fontsize", 23);
%icons = findobj(icons,'Type','line');
%icons = findobj(icons,'Marker','none','-xor');
%set(icons,'MarkerSize',15,'LineWidth',2);
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;

%ylabel("Beat contrast", "fontsize", 28)
%xlabel("Play speed (BPM)", "FontSize",28)
saveas(gca, "./summary/beat_contrast_4.png")
saveas(gca, "./summary/beat_contrast_4.fig")


%% plot the beat contrast of pos. 1 over pos. 2-4

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);
%temp = (bresall(:,:,1) - mean(bresall(:,:,2:4),3))./(bresall(:,:,1) + mean(bresall(:,:,2:4),3));
%errorbar((1:length(cands)),nanmean(temp,1),nanstd(temp,0,1)./sqrt(length(paths)),"ko-",'LineWidth',2)

temp = squeeze(mean(bresall(:,:,1)))./ squeeze(mean(mean(bresall(:,:,2:4),3)));
err = sqrt((nanstd(bresall(:,:,1),0,1)./squeeze(mean(mean(bresall(:,:,2:4),3)))).^2 ...
    + (nanstd(mean(bresall(:,:,2:4),3),0,1).*temp./squeeze(mean(mean(bresall(:,:,2:4),3)))).^2)/sqrt(size(bresall,1));
h1 = errorbar(1:length(cands), temp, err,"ko-",'LineWidth',3, "Markersize",15);
set(h1, 'markerfacecolor', get(h1, 'color'))
hold on

xticks(1:length(names))
xticklabels(cands + "%")
yticks(0.9:0.1:2)
ax = gca;
ax.YAxis.FontSize = 28;
ax.XAxis.FontSize = 28;
ylim([1.1,1.6])
%ylabel("pos.1 / pos.2-4", "fontsize", 28)
%xlabel("Play speed (BPM)", "FontSize",28)
xlim([0.5,length(names)+0.5])

saveas(gca, "./summary/beat_contrast_pos1.png")
saveas(gca, "./summary/beat_contrast_pos1.fig")

