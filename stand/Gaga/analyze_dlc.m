%% get sample beat
cd("/Users/TP/research/Tlab/stand/gaga")
dr = dir("./*LadyGaga*.wav");

[yda,Fs] = audioread(dr(1).folder + "/" + dr(1).name);

rng = [110, 131];
yda = yda(rng(1)*Fs:rng(2)*Fs,:);
xda = (1:length(yda))/Fs;

[yupper,ylower] = envelope(yda(:,1).*yda(:,1),1000,"peak");


%%
figure()
[pk,lc] = findpeaks(yupper,Fs,'MinPeakDistance',0.8);

plot(xda,yda(:,1).*yda(:,1))
hold on
plot(xda, yupper, "k")
plot(lc,pk,'kx', "linewidth",2)

pos = 6:13;
center =  6;
dif = mean(diff(lc(pos)))/2.015;
shift = 0.001;
beat = [];
for i = 0:60
    beat = [beat, lc(center)+shift + dif*i];
end

for i = 1:30
    beat = [beat, lc(center)+shift-dif*i];
end

beat = sort(beat);

beat(beat < -dif/2) = [];
beat(beat > (rng(2)- rng(1) + dif/2)) = [];
for b = 1:length(beat)
    plot([1,1]*beat(b),[0,0.5],"r");
end


%% the simultaneous plot of beat and DLC

csvpath = "./20210806_6_LadyGaga.csv";

csvname = dir(csvpath);
srate = 30;

dlc = readmatrix(csvname(1).folder + "/" +csvname(1).name);
data = dlc(rng(1)*srate:rng(2)*srate,:);


%nose = sqrt((data(:,2)-data(:,5)/2 - data(:,8)/2).^2 + (data(:,3)-data(:,6)/2 - data(:,9)/2).^2);
nose = sqrt((data(:,2)-mean(data(:,[5,8,11,14]),2)).^2 + (data(:,3)-mean(data(:,[6,9,12,15]),2)).^2);
%nose = diff(nose);

[pkd,lcd,width,prom] = findpeaks(nose,srate,'MinPeakProminence',2);
[pkd2,lcd2,width2,prom2] = findpeaks(-nose,srate,'MinPeakProminence',2);

thresh = 50;
pkd(prom > thresh) = [];
lcd(prom > thresh) = [];
prom(prom > thresh) = [];

pkd2(prom2 > thresh) = [];
lcd2(prom2 > thresh) = [];
prom2(prom2 > thresh) = [];

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.0, 0.6, 0.5]);
plot((data(1:end,1)-data(1,1))/srate,nose,"linewidth",2)
hold on
plot(lcd,pkd,'x','MarkerEdgeColor',[.8 0 .5],"linewidth",1.5)
plot(lcd2,-pkd2,'x','MarkerEdgeColor',[1 .5 0],"linewidth",1.5)

for i = 1:length(beat)
    plot([1,1]*beat(i),[0,300],"k");
end

xlim([0, rng(2) - rng(1)])
ylim([0,210])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;
ylabel(["nose displacement"],"Fontsize",26)
xlabel("time [s]","Fontsize",26)

%saveas(gcf, "./figure/gaga_findpeaks.png")


%%

figure()
datap = [];
for d = 1:length(lcd)
    [mn,I] = min(abs(beat - lcd(d)));
    datap = [datap, lcd(d)-beat(I)];
end
h1 = histogram(datap, -dif/2:dif/9:dif/2);
tmp = h1.Values;
datap2 = [];
for d = 1:length(lcd2)
    [mn,I2] = min(abs(beat - lcd2(d)));
    datap2 = [datap2, lcd2(d)-beat(I2)];
end
h2 = histogram(datap2, -dif/2:dif/9:dif/2);
tmp2 = [h2.Values;tmp];
hist = tmp2;
tmp2 = tmp2 / sum(sum(tmp2));

b = bar(tmp2', 1, 'stacked','FaceColor','flat');
b(1).CData = [1 .5 0];
b(2).CData = [.8 0 .5];
xlim([0,9]+0.5)
ylim([0,0.18])
xticks([0,9/4,9/2,9*3/4,9]+0.5)
xticklabels(["-π","-π/2","0","π/2","π"])

ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;

xlabel("Phase","Fontsize",26)
ylabel("Probability","Fontsize",26)

%saveas(gca, "./figure/gaga_histogram.png")

%%
%{
fig = figure();
datap = [];
for d = 1:length(lcd)
    [mn,I] = min(abs(beat - lcd(d)));
    datap = [datap, lcd(d)-beat(I)];
end
h1 = histogram(datap, -dif/2:dif/9:dif/2);
tmp = h1.Values/sum(h1.Values);
datap2 = [];
for d = 1:length(lcd2)
    [mn,I2] = min(abs(beat - lcd2(d)));
    datap2 = [datap2, lcd2(d)-beat(I2)];
end
h2 = histogram(datap2, -dif/2:dif/9:dif/2);
tmp2 = h2.Values/sum(h2.Values);

subplot(2,1,1)
b = bar(1:9,tmp,1,"stacked");
b.FaceColor = [0.8 0 .5];
xticks(1:2:9)
xticklabels(["-π","-π/2","0","π/2","π"])
ylim([0,0.25])
xlim([0,10])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;

subplot(2,1,2)
b2 = bar(1:9,tmp2,1,"stacked");
b2.FaceColor = [1 .5 0];
ylim([0,0.25])
xlim([0,10])
xticks(1:2:9)
xticklabels(["-π","-π/2","0","π/2","π"])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;
xlabel('Phase');

han=axes(fig,'visible','off'); 
ax = gca;
ax.YAxis.FontSize = 25;
han.YLabel.Visible='on';
ylabel(han,'Probability');
%saveas(gca, "./summary4/video_hist_human.png")

%}
%%

p = circ_rtest([datap,datap2]*2*pi/dif)  %0.0382

p = circ_rtest([datap]*2*pi/dif)  %0.2869

p = circ_rtest([datap2]*2*pi/dif)  %0.0979

