%% MUA plot of champion data


date = "210107";
path = "/Users/yoshiki/research/Tlab/data/IPC/IPC/";
dpath = "/Users/yoshiki/research/Tlab/data/IPC/IPC/raw/";
load(path + date + ".mat")

TC = ch_cortex;

trgindex = [9,7,2,4,5];  %counter = [1,2,2,2,1];

%%
mult = 32./[32,56,100,178,316];
labels = ["64","112","200","356","632"];

names = "datafile" + date + "_" + [32,56,100,178,316] + ".mat";
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.25, 1]);
for k = 1:length(mult)
    load(dpath + names(k))
    data = SpikeTiming;
    trg = trials{trgindex(k)}.click_or_silence;
    trg2 = trials{trgindex(k)}.trig_timings;
    disp([trg2(1)])
    st = 12500 + trg2(1);
    span = round(400/mult(k));
    rng = [st, st + span];
    resp = zeros(length(TC),max(trials{trgindex(k)}.trig_timings)+100000);
    ev = 10;
    for e = 1:length(TC)
        elenum = TC(e);
        for j = 1:ev
            resp(e,round(data{elenum}/30)+j) = resp(e,round(data{elenum}/30)+j) + 1;
        end
    end
    resp = resp(:,rng(1)+1:rng(2))*1000/ev;
    subplot(length(mult),1,k)
    [B,Is] = sort(mean(resp,2));
    mx = mean(max(resp,[],2));
    tmp = zeros(length(TC),span);
    for t = 1:length(TC)
        tmp(t,:) = resp(Is(t),:); 
        %plot(trg-1500,resp(Is(t),:)+t*mx,"-k")
        hold on
    end
    colormap(flipud(hot))
    imagesc(1:span,1:length(TC),tmp,[0,200])
    for i = 1:length(trg2)-1
        if trg2(i) > rng(1) && trg2(i) < rng(2) && trg(i) == 1
            plot([1,1]*(trg2(i)-st), [0,length(TC)+10], "-","color",[0.5,0.5,0.65],"linewidth",1.5)
            hold on
        end
    end
    rectangle('Position', [0,-5, 100, 10], "FaceColor", "g")
    xticks([])
    %xticklabels(13:0.5:(13+span/1000))
    yticks([1,length(TC)])
    yticklabels(["#1","#"+length(TC)])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    xlim([1,span])
    ylim([0, length(TC)])
    %h = colorbar("Fontsize",13);
    %ylabel(h, "MUA [s^-1]")
    %set(gca,'ytick',[])
    %set(gca,'yticklabel',[])
    title(labels(k)+ " [ms]","Fontsize",20)
end
%xlabel("time","Fontsize",20)

saveas(gca,path+"/figure/colormap_f.fig")
saveas(gca,path+"/figure/colormap_f.png")

