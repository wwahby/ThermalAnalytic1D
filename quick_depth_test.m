% dl = diag(U);
% dle = dl(end-1);

% strstr = 'foo (( bar (fie (fid)))) fee ( fie (foe (fum (foom)))) fur forg (fud (fow (foom)))';
% [depths, block_opens, block_closes] = find_block_depth(strstr);
[depths, block_opens, block_closes] = find_block_depth(detstr);

%%
figure(1)
clf
hold on
plot(depths)
% set(gca,'yscale','log')
% set(gca,'xscale','log')

%%
trigger_depth =  5;
all_inds = find(depths >= trigger_depth);
break_vec = [0, diff(all_inds)];
break_inds = find(break_vec > 1);
start_inds = [all_inds(1), all_inds(break_inds)];
stop_inds = [all_inds(break_inds-1), all_inds(end)]+1;

% start_inds = all_inds( mod(llv, 2) == 1);
% stop_inds = all_inds( mod(llv, 2) == 0);
% stop_inds = stop_inds + 1;
distances = stop_inds - start_inds;
dlen = length(distances);
dmax = max(distances);

func_cell = cell(1,length(start_inds));
str_cell = cell(1,2*length(start_inds)+1);
for ff = 1:length(start_inds)
    start_ind = start_inds(ff);
    stop_ind = stop_inds(ff);
    
%     func_str = ['@(a) ', detstr(start_ind:stop_ind) ];
%     func
    
    func_cell{ff} = str2func( ['@(a) ', detstr(start_ind:stop_ind) ]);
    
    str_cell{2*ff} = sprintf('( func_cell{%d}(a) )', ff);
    if ff > 1
        str_cell{2*ff-1} = detstr(stop_inds(ff-1)+1:start_ind-1);
    else
        str_cell{ff} = detstr(1:start_ind-1);
    end
end
str_cell{end} = detstr(stop_inds(end)+1:end);



newstr = strjoin(str_cell);
newfunc = str2func(newstr);

%%

figure(1)
clf
hold on
plot(depths)
% set(gca,'yscale','log')
% set(gca,'xscale','log')

% figure(2)
% clf
% hold on
% plot(depths(end-500:end))

%%
% figure(3)
% clf
% hold on
% plot(block_opens)
% plot(block_closes)
% xlim([length(block_opens)-500 length(block_opens)])
% 
% figure(4)
% clf
% hold on
% plot(block_opens)
% plot(block_closes)
% xlim([1 1000])

