function reg = regressERSP_Z(tfmap,RT, timewindow, times, freqs)

%%tfmap should be the single trial ERSP data freqsxtimesxtrials
%RT should be a row vector of values you want to regress, should have the same number of elements as number of trials 
%timewindow should have the latencies [min max] of the range to use
% times: timevector
% freqs: frequency vector

tfmaptmp = abs(tfmap);
baseind = find(times < 0);
timesind = find(times > timewindow(1) & times < timewindow(2));
tfmapsel = tfmaptmp(:,timesind,:);

baseline = repmat(mean(tfmaptmp(:,baseind,:),2),1,size(tfmapsel,2),1);

logtfmap = log10(tfmapsel)./log10(baseline);

Zmap=logtfmap;
ZRT = zscore(RT);

%% can use RTs as regressors 
%% single trial ERSPs

%% across trials - event related
for freq = 1:size(Zmap, 1)
    for time = 1:size(Zmap, 2)
        dep = squeeze(ZRT);
        ind = [ones([1 size(Zmap,3)]); squeeze(Zmap(freq,time,:))']; % cos and sin components of phase
        [b, bint, r, rint, stats] = regress(dep', ind'); clear bint r rint
        reg.R(freq,time) = stats(1);
        reg.F(freq,time) = stats(2);
        reg.p(freq,time) = stats(3);
        reg.errorvar(freq,time) = stats(4);
        reg.b(freq,time) = b(2); 
        clear b
        clear stats ind dep
    end
end

