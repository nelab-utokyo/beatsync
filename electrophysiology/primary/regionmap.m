%% creating corebelt segregation

paths = ["20200820","20200825","20200901","20200902","20200908","20200910","20200915","20200917","20201006","20201220","20201222","20210106","20210107","20210108_1","20210108_2"];

%load("electrode_position.mat")

lthresh = [11,11,12,11,10,10,11,10,10,11,12,12,11,12,11];
labels = [60,120,240,480];
ctrf = [];
resf = [];
resp = [];
rmvs = [];
p = 13;
pathf = "./"+paths(p) + "/figure/";
load("./"+paths(p)+"/shuf_param.mat")
load("./"+paths(p)+"/CF_latency.mat")

%TC = 1:96;
%figure()
%histogram(latency)

% available electrode map

%corebelt = zeros(96,1);
ampmap = zeros(10);
for i = 1:length(TC)
    ampmap(fig_mat_100(TC(i),1),fig_mat_100(TC(i),2)) = CF(TC(i))+1;
end

figure()
imAlpha=ones(size(ampmap));
imAlpha(isnan(ampmap))=0;
imagesc(ampmap,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);
for i = TC
    ps = fig_mat_100(i,:);
    %text(ps(2),ps(1),num2str(latency(i)))
    %if latency(i) < (lthresh(p) +1*(CF(i) < 10))
    if corebelt(i) == 1
        corebelt(i,1) = 1;
        rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"EdgeColor","r","Linewidth",1.5,"LineStyle","-")
        %if a1aaf(i) == 1
        %    rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"EdgeColor","g","Linewidth",1.5,"LineStyle","-")
        %else
        %    rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"EdgeColor","r","Linewidth",1.5,"LineStyle","-")
        %end
    end
end
tck = [1.6,2.0,2.5,3.2,4.0,5.0,6.4,8.0,10,13,16,20,25,32,40,50,57,64];
colorbar("Ticks",1:4:19,"Ticklabels",tck(2:4:19),"Fontsize",20)
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 18;
xlabel("Electrode position", "fontsize",28)
ylabel("Electrode position", "fontsize",28)
%}

%%
%
addpos = [7,2;6,2;5,2;4,2;3,3;1,4;8,7;1,6];
rmvpos = [8,10;8,9;7,9;6,9;3,9;3,8;6,7;3,6;10,9;10,8;2,1;2,2;];

for i = 1:size(addpos,1)
    num = find((fig_mat_100(:,1)==addpos(i,2))&(fig_mat_100(:,2)==addpos(i,1)));
    corebelt(num) = 1;
end

for i = 1:size(rmvpos,1)
    num = find((fig_mat_100(:,1)==rmvpos(i,2))&(fig_mat_100(:,2)==rmvpos(i,1)));
    corebelt(num) = 0;
end

figure()
imAlpha=ones(size(ampmap));
imAlpha(isnan(ampmap))=0;
imagesc(ampmap,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);

for i = TC
    ps = fig_mat_100(i,:);
    text(ps(2),ps(1),num2str(latency(i)))
    if corebelt(i) == 1
        rectangle("position", [ps(2)-0.5,ps(1)-0.5,1,1],"Linewidth",1.5,"LineStyle","--")
    end
end

ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 18;
xlabel("Electrode position", "fontsize",20)
ylabel("Electrode position", "fontsize",20)

%}
%%
pnt1 = [5,6];
pnt2 = [7,2];

a1aaf = zeros(96,1)-1;
for i = 1:10
    for j = 1:10
        num = find((fig_mat_100(:,1)==j)&(fig_mat_100(:,2)==i));
        if num > 96
            continue
        end
        if corebelt(num) == 1
            if i <= ((pnt2(1)-pnt1(1))*j + pnt2(2)*pnt1(1) - pnt2(1)*pnt1(2))/(pnt2(2)-pnt1(2))
            %if a1aaf(num) == 1
                a1aaf(num) = 1;
                text(i,j-0.2,"*")
            else
                a1aaf(num) = -1;
            end
        else
            a1aaf(num) = 0;
        end
    end
end
            
%%


save("./"+paths(p)+"/CF_latency.mat","corebelt","a1aaf","-append")


