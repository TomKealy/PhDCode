function copy_labels(g,h)
% copy_labels(g,h) --- copy labels from h to g

if ~is_labeled(h)
    error('Source graph (2nd argument) must be labeled')
end

ng = nv(g);
nh = nv(h);

label(g);
for v=1:min(ng,nh)
    label(g,v,get_label(h,v));
end