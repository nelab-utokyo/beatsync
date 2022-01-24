function [resp,reslfp] = data2matrix2(data,trg_,TC,ev)
    elenums = TC;
    span = 25;
    if ev <= 0
        resp = zeros(length(elenums),10,length(trg_));
        reslfp = zeros(length(elenums),10,length(trg_));

        for e = 1:length(elenums)
            elenum = TC(e);
            res = zeros(10,length(data.trg));
            res2 = zeros(10,length(data.trg));
            for i = 1:10
                res(i,data.("s"+i){elenum}+1) = res(i,data.("s"+i){elenum}+1) + 1;
                res2(i,:) = data.("d"+i)(elenum,:);
            end
            for j = 1:10
                for i = 1:length(trg_)
                    st = trg_(i) + 5;
                    resp(e,j,i) = mean(res(j,st:st+span-1));
                    %resp(e,j,i) = mean(res(j,st:st+span-1));
                    reslfp(e,j,i) = mean(res2(j,st:st+span-1));
                end
            end
        end
        resp = resp * 1000;
    elseif ev > 0
        resp = zeros(length(elenums),10,length(data.trg)+ev);
        reslfp = zeros(length(elenums),10,length(data.trg));

        for e = 1:length(elenums)
            elenum = TC(e);
            for i = 1:10
                for j = 1:ev
                    resp(e,i,data.("s"+i){elenum}+j) = resp(e,i,data.("s"+i){elenum}+j) + 1;
                end
                reslfp(e,i,:) = data.("d"+i)(elenum,:);
            end
        end
        resp(:,:,length(data.trg)+1:length(trg_)+ev) = [];
        resp = resp*1000/ev;
        if sum(abs(diff(trg_) - 1)) == 0
            resp = resp(:,:,trg_);
            reslfp = reslfp(:,:,trg_);
        end
    end
end

