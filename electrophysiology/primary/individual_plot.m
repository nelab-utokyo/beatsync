
%% raster plot of champion data 

champ = "20200917";
load("./"+champ+"/processed/exdata_p2f.mat")
load("./"+champ+"/shuf_param.mat")
load("./"+champ+"/buildmat/responsep.mat")
load("./"+champ+"/buildmat/responsef.mat")

%%
TC_ = TC;
TC_(rmv) = [];
mult = [1,2,4,8];
labels = ["60","120","240","480"];

names = "f" + [4,8,16,32]+"_15";
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.8]);
for k = 1:length(mult)

    data = df.(names(k));
    pkup = (21*mult(k)+1):(21*mult(k)+13);
    trg2 = data.trg2(pkup);
    trg = round(trg2(1)-100/mult(k)):(trg2(end)-1);
    disp([trg2(1),length(df.(names(k)).trg2)])
    [resp,resf] = data2matrix(data, trg, TC_, 1);
    subplot(4,1,k)
    resp = squeeze(mean(resp,2));
    [B,Is] = sort(mean(resp,2));
    mx = mean(max(resp,[],2));
    tmp = zeros(length(TC_),length(trg));
    for t = 1:length(TC_)
        tmp(t,:) = resp(Is(t),:); 
        %plot(trg-1500,resp(Is(t),:)+t*mx,"-k")
        hold on
    end
    colormap(flipud(hot))
    imagesc((1:length(trg))+trg(1)-1500,1:length(TC_),tmp,[0,200])
    for i = 1:length(trg2)-1
        plot([1,1]*(trg2(i)-1500), [0,length(TC_)+10], "-","color",[0.5,0.5,0.65],"linewidth",1.5)
        hold on
    end
    xticks([7000:500:11000])
    xticklabels(7:0.5:11)
    yticks([1,55])
    yticklabels(["#1","#55"])
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    xlim([trg(1)-1500,trg(end)-1500])
    ylim([0, length(TC_)])
    %h = colorbar("Fontsize",13);
    %ylabel(h, "MUA [s^-1]")
    %set(gca,'ytick',[])
    %set(gca,'yticklabel',[])
    title(labels(k)+ " BPM","Fontsize",25)
end
xlabel("time [s]","Fontsize",20)

saveas(gca,champ+"/figure/colormap_f.fig")
saveas(gca,champ+"/figure/colormap_f.png")



%% 
%%%%%%%%%%%%%%%%
% calculate musical on beat response
%%%%%%%%%%%%%%%%

paths = ["20200908","20200910","20200915","20200917","20201006","20210106","20210107"];

labels = [1,2,4,8]+"Hz";
splfp = 1;
ctrf = [];
cands = [75,100,200,300,400];

names = cands +"% (" + round((cands/100)*132) + "BPM)";

bresall = [];
nbresall = [];
for p = 1:length(paths)
    pathf = "./"+paths(p) + "/figure/";
    load("./"+paths(p)+"/shuf_param.mat")
    load("./"+paths(p)+"/buildmat/responsem.mat")
    load("./" + paths(p) + "/CF_latency.mat")
    load("./" + paths(p) + "/electrode_position.mat")
    TC(rmv) = [];
    bres_ = bres;
    bres_(:,rmv,:) = [];
    nbres_ = nbres;
    nbres_(:,rmv,:) = [];
    bresall_ = permute(bres_,[2,1,3]);
    nbresall_ = permute(nbres_,[2,1,3]);
    bresall = [bresall; bresall_];
    nbresall = [nbresall; nbresall_];
    
    contrast = (bresall_(:,:,5) - nbresall_(:,:,5))./(bresall_(:,:,5) + nbresall_(:,:,5));

    map = zeros(10,10,length(cands))*NaN;
    for j = 1:length(cands)
        counter = 0;
        for i = 1:length(TC)
            counter = counter + 1;
            map(fig_mat_100(TC(i),1),fig_mat_100(TC(i),2),j) = contrast(counter,j);
        end
    end
    %
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.59, 1]);
    poss = [[0.05 0.52 0.4 0.4];[0.5 0.52 0.4 0.4];[0.05 0.05 0.4 0.4];[0.5 0.05 0.4 0.4]];
    plt = [1,2,3,5];
    for i = 1:4
        subplot('Position',poss(i,:))
        tmp = squeeze(map(:,:,plt(i)));
        imAlpha=ones(size(tmp));
        imAlpha(isnan(tmp))=0;
        imagesc(tmp,'AlphaData',imAlpha,[0,0.65]);
        set(gca,'color',0*[1 1 1]);
        colorbar("Fontsize",18)
        title(names(plt(i)),"fontsize",30)
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        cnt = 0;
        for j = TC
            cnt = cnt + 1;
            if ismember(cnt,rmv)
                continue
            end
            ps = fig_mat_100(j,:);
            %if corebelt(j) == 1
            %    rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","r--")
            %end
            if a1aaf(j) == 1.5
                %text(ps(2)-0.1,ps(1),"x")
                rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","-","EdgeColor",[0.4,1,0.4])
            elseif abs(a1aaf(j)) == 1
                %text(ps(2)-0.1,ps(1),"o")
                rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","-","EdgeColor",[1,0.3,0.3])
            end
        end  
    end
    saveas(gca, pathf + "music_contrast_map.png")
    %}
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.45, 0.6]);
    h = boxplot(contrast);
    set(h,'LineWidth',2)
    hold on
    for i = 1:size(contrast,2)
        scatter(rand(size(contrast,1),1)*0.08+i-0.04, contrast(:,i),130,".","markeredgecolor",[0.1,0.5,0])
    end
    xticks(1:length(names))
    xticklabels(cands + "%")
    yticks(-0.4:0.2:1)
    %ylim([0.3,0.5])
    xlim([0.5,length(names)+0.5])
    ax = gca;
    ax.YAxis.FontSize = 28;
    ax.XAxis.FontSize = 28;
    %ylabel("Beat contrast", "FontSize",28)
    %xlabel("Play speed(BPM)", "FontSize",28)
    saveas(gca, pathf + "music_contrast.png")
end


%% position map of beat contrast (champion data) 

champ = "20210107";
load("./"+champ+"/shuf_param.mat")
load("./"+champ+"/buildmat/responsem.mat")
load("./"+champ+"/CF_latency.mat")
load("./"+champ+"/electrode_position.mat")

%%

contrast = (bresall_(:,:,5) - nbresall_(:,:,5))./(bresall_(:,:,5) + nbresall_(:,:,5));

contrast(rmv,:) = [];
TC(rmv) = [];
map = zeros(10,10,size(contrast,2))*NaN;
for j = 1:size(contrast,2)
    for i = 1:length(TC)
        %if sum(rmv == i) > 0
        %    continue
        %end
        map(fig_mat_100(TC(i),1),fig_mat_100(TC(i),2),j) = contrast(i,j);
    end
end

names = cands +"%(" + round((cands/100)*132) + ")";

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.5]);
%poss = [[0.05 0.52 0.4 0.4];[0.5 0.52 0.4 0.4];[0.05 0.05 0.4 0.4];[0.5 0.05 0.4 0.4]];
counter = 0;
for i = [1,3,5]
    %subplot('Position',poss(i,:))
    counter = counter + 1;
    subplot(1,3,counter)
    tmp = squeeze(map(:,:,i));
    imAlpha=ones(size(tmp));
    imAlpha(isnan(tmp))=0;
    imagesc(tmp,'AlphaData',imAlpha,[0.,0.75]);
    set(gca,'color',0*[1 1 1]);
    colorbar('XTick', 0:0.25:0.75, "Fontsize",18)
    title(names(i),"fontsize",30)
    set(gca,'YTickLabel',[])
    set(gca,'XTickLabel',[])
    cnt = 0;
    for j = TC
        cnt = cnt + 1;
        if ismember(cnt,rmv)
            continue
        end
        ps = fig_mat_100(j,:);
        %if corebelt(j) == 1
        %    rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","r--")
        %end
        if a1aaf(j) == 1
            %text(ps(2)-0.1,ps(1),"x")
            rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","-","EdgeColor",[0.4,1,0.4])
        elseif a1aaf(j) == -1
            %text(ps(2)-0.1,ps(1),"o")
            rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","-","EdgeColor",[1,0.3,0.3])
        end
    end  
end

saveas(gca,"./"+champ + "/figure/mucic_contrast_map.fig")
saveas(gca,"./"+champ + "/figure/music_contrast_map.png")  
