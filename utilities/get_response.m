function Sr = get_response(S, ti, smooth)
    if smooth
        DFF = S.sDFF;
    else
        DFF = S.DFF;
    end
    
    K = S.K;
    
    baseline_STD = NaN(K, ti.num_stim);
    baseline_MED = NaN(K, ti.num_stim);
    peak_response = NaN(K, ti.num_stim);
    peak_ind = NaN(K, ti.num_stim);

    for it = 1:ti.num_stim
        sdev=std(DFF(:,ti.baseline_start(it):ti.baseline_end(it)),1,2);
        baseline_STD(:,it)=sdev;

        med=median(DFF(:,ti.baseline_start(it):ti.baseline_end(it)),2);
        baseline_MED(:,it) = med;

        %peak_response(:,it)=max(DFF(:,ti.peakresp_start(it):ti.peakresp_end(it)),[],2);
        [peak_response(:,it), peak_ind(:,it)]=max(DFF(:,ti.peakresp_start(it):ti.peakresp_end(it)),[],2);
    end
    peakAmp = peak_response - baseline_MED;
    response = (peakAmp>(baseline_MED+3*baseline_STD)) & (peakAmp>1);
    peak_ind = peak_ind + ti.peakresp_start;
    
    Sr = S;
    Sr.peak_ind=peak_ind;
    Sr.baseline_STD = baseline_STD;
    Sr.baseline_MED = baseline_MED;
    Sr.peak_response = peak_response;
    Sr.peakAmp = peakAmp;
    Sr.response = response;
    
end