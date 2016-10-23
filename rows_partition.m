m = 100;
n = 200;

m_p = m/50;

A = randn(m,n);
offset = 1;
for i=1:50
   rows = [offset:offset+(m_p-1)]
   Ap = A(rows, :);
   offset = offset + m_p;
end;