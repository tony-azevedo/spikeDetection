function all_filtered_data = filterDataWithSpikes(vars)

if ~isempty(vars.spikeTemplate)
    vars.spikeTemplateWidth = length(vars.spikeTemplate);
end
filts1 = vars.hp_cutoff/(vars.fs/2);
[x,y] = butter(3,filts1,'high');%%bandpass filter between 50 and 200 Hz
filtered_data_high = filter(x, y, vars.unfiltered_data-vars.unfiltered_data(1));

filts2 = vars.lp_cutoff/(vars.fs/2);
[x2,y2] = butter(3,filts2,'low');%%bandpass filter between 50 and 200 Hz
filtered_data = filter(x2, y2, filtered_data_high);
filtered_data = filtered_data(:)';

if vars.diff == 0
    diff_filt = filtered_data;
elseif vars.diff == 1
    diff_filt = [0 diff(filtered_data)];
    diff_filt(1:100) = 0;
elseif vars.diff == 2
    diff_filt = [0 0 diff(diff(filtered_data))'];
    diff_filt(1:100) = 0;
end

all_filtered_data = diff_filt;

all_filtered_data = vars.polarity*all_filtered_data;


