% quick_eval_test

aa_str = { '@(a) 1+a', '@(a) 1-a'};
aa_func = { @(a) 1+a, @(a) 1-a};

funcstr = '@(a) 10*aa_func{1}(a)';
ff = eval(funcstr);

ff(1)
ff(2)

gg = str2func(funcstr);
gg(1)
gg(2)