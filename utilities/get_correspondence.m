function S = get_correspondence(binC)
    [X,Y] = size(binC);
    X=1:X;
    Y=1:Y;
    
    [r,c] = find(binC);
    corresp_list = [r c];

    unique_labels = unique(r);
    counts = sum(binC,2);
    value_counts = [(1:size(binC,1))', counts];
    one_to_one = find(counts==1);
    repeats = find(counts>1);
    
    idx = ismember(r, one_to_one);
    uniquely_assigned = corresp_list(idx,:);
    assigned_X = uniquely_assigned(:,1);
    assigned_Y = uniquely_assigned(:,2);
    unassigned_X = setdiff(X,assigned_X);
    unassigned_Y = setdiff(Y,assigned_Y);
    
    S.binC = binC;
    S.corresp_list = corresp_list;
    S.unique_labels = unique_labels;
    S.value_counts = value_counts;
    S.one_to_one = one_to_one;
    S.repeats = repeats;
    S.uniquely_assigned = uniquely_assigned;
    S.assigned_X = assigned_X;
    S.assigned_Y = assigned_Y;
    S.unassigned_X = unassigned_X;
    S.unassigned_Y = unassigned_Y;
    
end