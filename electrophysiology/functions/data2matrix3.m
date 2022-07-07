function [resp,reslfp,ssum] = data2matrix3(data,trg,TC,ev)
    elenums = TC;
    span = 25;
    if ev <= 0
        trg_ = [trg, trg(end)+median(diff(trg))*3];
        %disp(median(diff(trg))*3)
        resp = zeros(length(elenums),10,length(trg));
        reslfp = zeros(length(elenums),10,length(trg));
        ssum = zeros(length(elenums)+1,1);
        for e = 1:length(elenums)
            elenum = TC(e);
            res = zeros(10,length(data.trg));
            res2 = zeros(10,length(data.trg));
            for i = 1:10
                res(i,data.("s"+i){elenum}+1) = res(i,data.("s"+i){elenum}+1) + 1;
                res2(i,:) = data.("d"+i)(elenum,:);
            end
            for j = 1:10
                for i = 1:length(trg)
                    st = trg_(i) + 5;
                    en = st + span - 1;
                    resp(e,j,i) = mean(res(j,st:en)) - mean(res(j,trg_(i):st -1));
                    reslfp(e,j,i) = mean(res2(j,st:en)) - mean(res2(j,trg_(i):st-1));
                end
            end
            ssum(e,1) = sum(mean(res(:,trg(1):trg(end)),1));
        end
        ssum(end,1) = trg(end) - trg(1) + 1;
        resp = resp * 1000;
    elseif ev > 0
        resp = zeros(length(elenums),10,length(data.trg)+ev);
        reslfp = zeros(length(elenums),10,length(data.trg));
        ssum = zeros(length(elenums),1);
        for e = 1:length(elenums)
            elenum = TC(e);
            for i = 1:10
                for j = 1:ev
                    resp(e,i,data.("s"+i){elenum}+j) = resp(e,i,data.("s"+i){elenum}+j) + 1;
                end
                reslfp(e,i,:) = data.("d"+i)(elenum,:);
            end
            ssum(e,1) = sum(mean(resp(:,trg(1):trg(end)),1));
        end
        resp(:,:,length(data.trg)+1:length(data.trg)+ev) = [];
        resp = resp*1000/ev;
        if sum(abs(diff(trg) - 1)) == 0
            resp = resp(:,:,trg);
            reslfp = reslfp(:,:,trg);
        end
    end
end

