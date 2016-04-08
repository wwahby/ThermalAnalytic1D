% diagU = diag(U);
% diagL = diag(L);
% 
% detA = det(Pprod(diag(L))*prod(diag(U));
% 
% fprintf('Converting syms...')
% ta = cputime;
% detstr = char(detA);
% tb = cputime;
% fprintf('\t%.3g s\n', tb- ta);

%%
fprintf('Counting parens...')
ta = cputime;
aaa = strfind(detstr, '(');
bbb = strfind(detstr, ')');
tb = cputime;
fprintf('\t%.3g s\n', tb- ta);

fprintf('Making functions...')
ta = cputime;
% block_opens = @(ii) find(aaa < ii, 1, 'last');
% block_closes = @(ii) find( bbb < ii, 1, 'last');

block_depth = @(ii) block_opens(ii) - block_closes(ii);
tb = cputime;
fprintf('\t%.3g s\n', tb- ta);

%%

% depths = zeros(1,max(bbb));
% depths(aaa(1):bbb(1)-1) = 1;

fprintf('Making vectors...')
ta = cputime;
block_opens = zeros(1,max(bbb));
block_closes = zeros(1,max(bbb));
for iind = 1:length(aaa)-1
    if mod(iind,1e4) == 0
        fprintf('%d / %d\n', iind, max(bbb) );
    end
    block_opens( aaa(iind):aaa(iind+1)-1 ) = iind;
    block_closes( bbb(iind):bbb(iind+1)-1 ) = iind;
end
block_opens(aaa(iind+1):end) = iind+1;
block_closes(bbb(iind+1):end) = iind+1;
tb = cputime;
fprintf('\t%.3g s\n', tb- ta);

%%
depths = block_opens - block_closes;


figure(1)
clf
hold on
plot(depths)
% set(gca,'yscale','log')