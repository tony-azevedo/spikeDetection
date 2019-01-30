function [vars,spikeWaveform,spikeWaveform_] = likelyInflectionPoint(vars,spikeWaveforms,targetSpikeDist)

% normalize and find a peak of the second derivative
idx_i = round(vars.spikeTemplateWidth/6);
idx_f = round(vars.spikeTemplateWidth/24);
idx_m = round(vars.spikeTemplateWidth*2/3);

window = -floor(vars.spikeTemplateWidth/2): floor(vars.spikeTemplateWidth/2);
spikewindow = window-floor(vars.spikeTemplateWidth/2);
smthwnd = (vars.fs/2000+1:length(spikewindow)-vars.fs/2000);

% Find the best estimate of the spike shape and it's 2nd derivative
for sp = 1:size(spikeWaveforms,2)
    spikeWaveforms(:,sp) = spikeWaveforms(:,sp)-min(spikeWaveforms(:,sp));
    spikeWaveforms(:,sp) = spikeWaveforms(:,sp)/max(spikeWaveforms(:,sp));
end

goodspikes = targetSpikeDist<quantile(targetSpikeDist,.25);
if sum(goodspikes)<4
    goodspikes(:) = 1;
end
spikeWaveform = nanmean(spikeWaveforms(:,goodspikes),2);
spikeWaveform = smooth(spikeWaveform-spikeWaveform(1),vars.fs/2000);
spikeWaveform_ = smoothAndDifferentiate(spikeWaveform,vars.fs/2000);
spikeWaveform_ = spikeWaveform_-spikeWaveform_(smthwnd(1));
spikeWaveform_ = (spikeWaveform_-min(spikeWaveform_(idx_i:end-idx_f)))/diff([min(spikeWaveform_(idx_i:end-idx_f)) max(spikeWaveform_(idx_i:end-idx_f))]);

% find the peak of the second derivative that is not the end point
[pks,inflPntPeak_ave] = findpeaks(spikeWaveform_(idx_i+1:end-idx_f),'MinPeakProminence',0.014);
inflPntPeak_ave = inflPntPeak_ave+idx_i;
pks = pks(abs(inflPntPeak_ave-idx_m)==min(abs(inflPntPeak_ave-idx_m)));
inflPntPeak_ave = inflPntPeak_ave(abs(inflPntPeak_ave-idx_m)==min(abs(inflPntPeak_ave-idx_m)));
if numel(inflPntPeak_ave) > 1 
    % if by chance the peaks are equidistant from the estimated mid point
    inflPntPeak_ave = inflPntPeak_ave(pks==max(pks));
end
vars.likelyiflpntpeak = inflPntPeak_ave; 