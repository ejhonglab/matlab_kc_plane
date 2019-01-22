function f = plotCorrMat(response,ti,range, tt)
% PLOTCORRMAT   plots correlation matrices given an array of responses (# 
%               components x stimuli, in order that they were presented,
%               then rearranged so they are grouped by odor)

% Inputs
%   response:   [# component x # stim]  array, each column is a response vector
%   ti:         "trigger info", timestamp information + stimulus + frame numbers
%   range:       range of colormap to display
%   tt:          optional, title figure

% Outputs
%   f:          figure handle
%

if ~exist('range', 'var')
   range = [-1 1];
end

if ~exist('tt', 'var')
    tt = 'Response corr';
end

K = size(response,1);               % # of components
[out,idx] = sort(ti.si);            % idx = indices used to sort columns by odor


%shuffle responses 
response_shuffled = NaN(size(response));
for i = 1:size(response,2)
    response_shuffled(:,i) = response(randperm(K),i);        
end
response_sorted = response(:,idx);
response_sorted_shuffled = response_shuffled(:,idx);

%% calculate correlation coefficients between response vectors
R = corrcoef(response);             % presentation order
Ro = corrcoef(response_sorted);     % grouped by odor
R_shuffled = corrcoef(response_shuffled);   % shuffle control, presentation order
Ro_shuffled = corrcoef(response_sorted_shuffled); % shuffle control,grouped by odor

%% plot correlation matrices
f = figure('Name', 'Response corr. matrices');
dmap=flipud(cbrewer('div','RdBu', 100));
colormap(dmap)

% calculate locations to plot axis labels
block = ones(1,ti.num_stim);
for i =1: ti.num_stim    
    temp = (ti.stim_ict(i)>ti.block_ict & ti.stim_ict(i)<ti.block_fct)...
        & (ti.stim_fct(i)>ti.block_ict & ti.stim_fct(i)<ti.block_fct);
    temp = find(temp);
    block(i) = temp;
end
li = find(diff(block)==1)+.5;
block_centers = [];
for i = 1:numel(unique(block))
    block_centers(i) = mean(find(block==i));
end

% calculate locations to draw lines separating odors
li = find(diff(out)==1)+.5;
[U,ia,ic] = unique(out);
a_counts = accumarray(ic,1);
oi = cumsum(a_counts)';
centers = [];
for i = 1:numel(unique(ic))
    centers(i) = mean(find(ic==i));
end

subplot(2,2,1) % top left corner, presentation order
    h=imagesc(R, range);
    h.AlphaData = ones(size(h.CData)); 
    h.AlphaData(isnan(h.CData)) = 0;
    axis image; colorbar;
    %set(gca, 'Color', 'k')    
    title('Presentation order');
    set(gca, 'XTick', block_centers, 'XTickLabel', 1:ti.num_trials, 'XAxisLocation', 'top',...
       'YTick', block_centers, 'YTickLabel', 1:ti.num_trials, 'TickLabelInterpreter', 'none');
    vline(li, 'k-');
    hline(li, 'k-');
subplot(2,2,2) % top right corner, shuffle control in presentation order 
    h=imagesc(R_shuffled, range);
    h.AlphaData = ones(size(h.CData)); 
    h.AlphaData(isnan(h.CData)) = 0;
    axis image; colorbar;
    %set(gca, 'Color', 'k')
    title('Presentation order, shuffle control');
     set(gca, 'XTick', block_centers, 'XTickLabel', 1:ti.num_trials, 'XAxisLocation', 'top',...
       'YTick', block_centers, 'YTickLabel', 1:ti.num_trials, 'TickLabelInterpreter', 'none');
    vline(li, 'k-');
    hline(li, 'k-');
    vline(li, 'k-');
    hline(li, 'k-');

subplot(2,2,3) % bottom left, grouped by odor
    h=imagesc(Ro, range);
    h.AlphaData = ones(size(h.CData)); 
    h.AlphaData(isnan(h.CData)) = 0;
    axis image; colorbar;
    %set(gca, 'Color', 'k')
    title('Sorted by odor');
    set(gca, 'XTick', centers, 'XTickLabel', ti.pin_odors, 'XAxisLocation', 'top',...
       'YTick', centers, 'YTickLabel', ti.pin_odors, 'TickLabelInterpreter', 'none');
    ax=gca;
    ax.YAxis.FontSize=8;
    ax.XAxis.FontSize=8;
    xtickangle(45);
    vline(li, 'k-');
    hline(li, 'k-');

subplot(2,2,4) % bottom right, grouped by odor
    h=imagesc(Ro_shuffled, range);
    h.AlphaData = ones(size(h.CData)); 
    h.AlphaData(isnan(h.CData)) = 0;
    axis image; colorbar;
    %set(gca, 'Color', 'k')
    title(['Sorted by odor, shuffle control']);
    set(gca, 'XTick', centers, 'XTickLabel', ti.pin_odors, 'XAxisLocation', 'top',...
       'YTick', centers, 'YTickLabel', ti.pin_odors, 'TickLabelInterpreter', 'none');
    ax=gca;
    ax.YAxis.FontSize=8;
    ax.XAxis.FontSize=8;
    xtickangle(45);
    vline(li, 'k-');
    hline(li, 'k-');

tt = suptitle(tt);
tt.Interpreter = 'none';
tt.FontSize = 10;

end

