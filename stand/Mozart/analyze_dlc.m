%% get sample beat

cd("/Users/yoshiki/research/Tlab/stand/20210806_mozart_3")
dr = dir("./*Mozart*.wav");

[yda,Fs] = audioread(dr(1).folder + "/" + dr(1).name);

rng = [90, 120];
yda = yda(rng(1)*Fs:rng(2)*Fs,:);
xda = (1:length(yda))/Fs;

[yupper,ylower] = envelope(yda(:,1).*yda(:,1),1000,"peak");


%%
figure()
[pk,lc] = findpeaks(yupper,Fs,'MinPeakDistance',0.8);

plot(xda,yda(:,1).*yda(:,1))
hold on

%plot(lc,pk,'x')

dif = mean(diff(lc(end-6:end-4)))/2.005;

beat = [];
for i = 0:13
    plot([1,1]*(lc(end-5)+0.02+dif*i),[0,0.5],"r");
    beat = [beat, lc(end-5)+dif*i];
end

for i = 1:51
    plot([1,1]*(lc(end-5)+0.02-dif*i),[0,0.5],"r");
    beat = [beat, lc(end-5)-dif*i];
end

beat = sort(beat);


%% the simultaneous plot of beat and DLC

csvpath = "./*csv";

csvname = dir(csvpath);

dlc = readmatrix(csvname(1).folder + "/" +csvname(1).name);

srate = 30;

data = dlc(rng(1)*srate:rng(2)*srate,:);
%nose = sqrt((data(:,2)-data(:,5)/2 - data(:,8)/2).^2 + (data(:,3)-data(:,6)/2 - data(:,9)/2).^2);
nose = sqrt((data(:,2)-mean(data(:,[5,8,11,14]),2)).^2 + (data(:,3)-mean(data(:,[6,9,12,15]),2)).^2);

%[pkd,lcd] = findpeaks(nose,srate,'MinPeakProminence',2.5);
%[pkd2,lcd2] = findpeaks(-nose,srate,'MinPeakProminence',2.5);
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
plot(lcd,pkd,'x','MarkerEdgeColor',[.8 0 .5], "linewidth",1.5)
plot(lcd2,-pkd2,'x','MarkerEdgeColor',[1 .5 0],"linewidth",1.5)

for i = 1:length(beat)
    plot([1,1]*beat(i),[0,210],"k");
end

ylim([0,210])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;
ylabel(["nose displacement"],"Fontsize",26)
xlabel("time [s]","Fontsize",26)

saveas(gcf, "./figure/Mozart_findpeaks.png")


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
tmp2 = hist / sum(sum(tmp2));

b = bar(tmp2', 1, 'stacked','FaceColor','flat');
b(1).CData = [1 .5 0];
b(2).CData = [.8 0 .5];
xlim([0,9]+0.5)
xticks([0,9/4,9/2,9*3/4,9]+0.5)
xticklabels(["-π","-π/2","0","π/2","π"])
ax = gca;
ax.YAxis.FontSize = 23;
ax.XAxis.FontSize = 23;
ylim([0,0.2])
xlabel("Phase","Fontsize",26)
ylabel("Probability","Fontsize",26)

saveas(gca, "./figure/Mozart_histogram.png")


%%
figure()
polarscatter([datap,datap2]*2*pi/dif,[prom;prom2]',"filled")
p = circ_rtest([datap,datap2]*2*pi/dif,[prom;prom2]')  %6.55e-58


%% rayleight's test

p = circ_rtest([datap,datap2]*2*pi/dif)  %0.0298


%%

alpha = [datap,datap2]*2*pi/dif;
alpha = circ_rad2ang(alpha);
n = length(alpha);
alpha = sort(alpha);

U = 0;
lambda = 360/n;
for j = 1:n-1
    ti = alpha(j+1) - alpha(j);
    U = U + abs(ti - lambda);
end
U = (1/2)*U;


%% kai squared test

allhist = sum(hist,1);
m = sum(allhist)/9;

kai = sum((allhist - m).^2)/m;
disp(kai)


