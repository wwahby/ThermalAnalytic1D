function [newfunc, newstr, func_cell] = condense_functions(old_str, trigger_depth)
[depths, ~, ~] = find_block_depth(old_str);
all_inds = find(depths >= trigger_depth);
if (~isempty(all_inds))
    break_vec = [0, diff(all_inds)];
    break_inds = find(break_vec > 1);
    start_inds = [all_inds(1), all_inds(break_inds)];
    stop_inds = [all_inds(break_inds-1), all_inds(end)]+1;

    func_cell = cell(1,length(start_inds));
    str_cell = cell(1,2*length(start_inds)+1);
    for ff = 1:length(start_inds)
        start_ind = start_inds(ff);
        stop_ind = stop_inds(ff);
        
        func_cell{ff} = eval( ['@(a) ', old_str(start_ind:stop_ind) ]); % can't use str2func for this after 2015a

        str_cell{2*ff} = sprintf('( func_cell{%d}(a) )', ff);
        if ff > 1
            str_cell{2*ff-1} = old_str(stop_inds(ff-1)+1:start_ind-1);
        else
            str_cell{ff} = old_str(1:start_ind-1);
        end
    end
    str_cell{end} = old_str(stop_inds(end)+1:end);

    newstr = strjoin(str_cell);
    newfunc = eval(newstr); % can't use str2func for this after 2015a
else
    newstr = old_str;
    newfunc = str2func(newstr);
    func_cell = {};
end

