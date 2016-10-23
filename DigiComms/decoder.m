function y_dec=decoder(y,code_table)
for m=1:size(code_table,1)
   if y==code_table(m,:), y_dec=m; break; end
end
