function th=principal_frequency(th,thl,thu)
if nargin<3, thu=pi; end
if nargin<2, thl=-pi; end
while sum(th<=thl)>0, ind=find(th<=thl); th(ind)=th(ind)+2*pi; end
while sum(th>thu)>0, ind=find(th>thu); th(ind)=th(ind)-2*pi; end 