function p = find_path(g,u,v)
% find_path(g,u,v) --- find a shortest path from u to v

n = nv(g);

if (u<1) || (u>n) || (v<1) || (v>n)
    p = [];
    return
end

if u==v
    p = u;
    return
end

q_init(n+1); 
track = zeros(1,n);
track(v) = v;

q_push(v); 

while(q_size > 0)
    t = q_pop_front;
    if t==u
        break
    end
    push_list = neighbors(g,t);
    for s = push_list
        if (track(s) == 0)
            track(s) = t;
            q_push(s);
        end
    end
end

if track(u) == 0
    p = [];
    return
end

p = [];
last = u;
while (last ~= v)
    p = [p,last];
    last = track(last);
end
p = [p,v];


q_init(1);
