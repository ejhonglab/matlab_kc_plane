function response = compute_response_vectors(DFF,ti,responseOptions)
%UNTITLED2 Summary of this function goes here
% Inputs:
%   DFF                   dF/F traces, one row per neuron/component
%   ti                  trigger info, need ti.num_stim, ti.baseline_start/end,
%                          ti.peakresp_start/end
%   responseOptions     optional
%                          .std_thr, threshold for stdev around baseline
%                          .thr, threshold for change in dF/F from baseline
% 
% 
% 
% 


% Outputs:
%      response
%       .baseline_STD
%       .baseline_MED
%       .peak_response
%       .peakAmp
%       .bin_response
%       .peak_ind
% 
% 


if ~exist('responseOptions', 'var')
    responseOptions.std_thr = 3;
    responseOptions.thr = 1;
end

K = size(DFF,1);
    
baseline_STD = NaN(K, ti.num_stim);
baseline_MED = NaN(K, ti.num_stim);
peak_response = NaN(K, ti.num_stim);
peak_ind = NaN(K, ti.num_stim);

for it = 1:ti.num_stim
    sdev=std(DFF(:,ti.baseline_start(it):ti.baseline_end(it)),1,2);
    baseline_STD(:,it)=sdev;

    med=median(DFF(:,ti.baseline_start(it):ti.baseline_end(it)),2);
    baseline_MED(:,it) = med;

    [peak_response(:,it), peak_ind(:,it)]=max(DFF(:,ti.peakresp_start(it):ti.peakresp_end(it)),[],2);
end
peakAmp = peak_response - baseline_MED;
bin_response = (peakAmp>(baseline_MED+responseOptions.std_thr*baseline_STD)) & (peakAmp>responseOptions.thr);
peak_ind = peak_ind + ti.peakresp_start;

response.peak_ind=peak_ind;
response.baseline_STD = baseline_STD;
response.baseline_MED = baseline_MED;
response.peak_response = peak_response;
response.peakAmp = peakAmp;
response.bin_response = bin_response;
response.options = responseOptions;
end

