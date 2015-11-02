%% VERSIONE 20 11 2011
function s_out = repeat_n(s_in,n)

s_out = s_in;
for j = 1:(n-1)
    s_out = [s_out s_in(2:end)];
end
