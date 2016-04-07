function str_c = convert_strfunc_for_func(astr)
str_a = strrep( astr, '*', '.*');
str_b = strrep(str_a, '^', '.^');
str_c = strrep(str_b, '/', './');