function [ft] = plot_S_traces(S,ti,cids,response,tt)
%PLOT_TRACES Summary of this function goes here
%   response:       response struct (response.bin_response used to plot)
%   tt:             plot title (optional)
%   

if ~exist('tt', 'var')
    tt = 'dF/F traces';
end

% colors
c1 = [0.9294 0.6902 0.1294];
c2 = [0.8889    0.9444    0.4000];
c3 = [0.9290 0.6940 0.1250];

% plot traces
ft = figure;
Hmw = iosr.figures.multiwaveplot(ti.frame_times, [], S.sDFF(cids,:),...
    'mode', 'plot',...
    'reverseY', true,...
    'gain', 3);

% draw patches to show stimulus delivery
yl = [ylim() fliplr(ylim())];
y = repmat(yl, ti.num_stim, 1)';
x =  [ti.stim_ict ti.stim_ict ti.stim_fct ti.stim_fct]';
patch(x,y,c2, 'FaceAlpha', .4, 'LineStyle', 'none');    

% label stimuli
hold on
set(gca, 'XTick', ti.frame_times(ti.stim_on),...
    'XTickLabel', ti.pin_odors(ti.si),...
    'XAxisLocation', 'top');
ax=gca;
ax.YAxis.FontSize=8;
ax.XAxis.FontSize=8;
xtickangle(45);

title(tt, 'Interpreter', 'None')
    

% add peaks if 'response' was an input variable
if exist('response', 'var')
    [r,c] = find(response.bin_response(cids,:));
    x= nan(size(r));
    y = nan(size(c));
    for i = 1:numel(r)
        ri = r(i);
        ci = response.peak_ind(cids(ri),c(i));
        x(i) = Hmw(ri).XData(ci);
        y(i) = Hmw(ri).YData(ci);
    end
    scatter(x,y-0.15, 10,'ro', 'filled');
end

for i = 1:ti.num_trials-1
    y = [ylim() fliplr(ylim())];
    x =  [ti.frame_times(ti.trial_end(i)) ti.frame_times(ti.trial_end(i)),...
    ti.frame_times(ti.trial_start(i+1)) ti.frame_times(ti.trial_start(i+1))]';    
    patch(x,y,'w','LineStyle', 'none');    

end

end

