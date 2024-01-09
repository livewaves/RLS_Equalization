function en = env (x, L);

x_length  = length (x);
for i = L:(x_length-L);
    en (i) = max (x (i:(i+L)));
end;
en (1:(L-1)) = max (x (1:L));