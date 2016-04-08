function [depths, block_opens, block_closes] = find_block_depth(funcstr)
%%
fprintf('Counting parens...')
ta = cputime;
aaa = strfind(funcstr, '(');
bbb = strfind(funcstr, ')');
tb = cputime;
fprintf('\t%.3g s\n', tb- ta);

%%
fprintf('Making vectors...')
ta = cputime;
block_opens = zeros(1,max(bbb));
block_closes = zeros(1,max(bbb));
for iind = 1:length(aaa)-1
    block_opens( aaa(iind):aaa(iind+1)-1 ) = iind;
    block_closes( bbb(iind):bbb(iind+1)-1 ) = iind;
end
block_opens(aaa(iind+1):end) = iind+1;
block_closes(bbb(iind+1):end) = iind+1;
tb = cputime;
fprintf('\t%.3g s\n', tb- ta);

%%
depths = block_opens - block_closes;

% 
% figure(1)
% clf
% hold on
% plot(depths)
% % set(gca,'yscale','log')
% set(gca,'xscale','log')