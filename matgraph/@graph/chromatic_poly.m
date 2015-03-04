function out = chromatic_poly(g)
% chrompoly(g) --- find the chromatic polynomial of g
% Warning: This algorithm is slow and unusable except for small graphs.
% Author: James Preen

out = chrompoly(double(matrix(g)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is James Preen's .m file very lightly edited %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = chrompoly(A)

out=[];
n=size(A,1);
e=sum(sum(A))/2;

if e==0,
    out(n+1)=0;out(1)=1;
	%  disp(['that is k bar ' num2str(n)]);
else
	somenew=1;
	nb=ceil(rand(1)*n);
	while somenew,
		nb2=rem(find(A(nb',:)')-1,n)+1;
		nb2=nb2';
		nb3=unique([nb nb2]);
		if size(nb3,2)==size(nb,2), somenew=0;end;
		nb=nb3;
	end
	if size(nb,2)<n,
		%       disp(['disconnected set: ' num2str([nb])]);
		nb2=setdiff(1:n,nb);
		A1=A;A1(:,nb2)=[];A1(nb2,:)=[];
		A2=A;A2(:,nb)=[];A2(nb,:)=[];
		out=conv(chrompoly(A1),chrompoly(A2));
	else
		valseq=sum(A);
		% v1 = min(find(valseq==max(valseq)));
		v1 = find(valseq==max(valseq),1);
		if valseq(v1)==n-1,
			%   disp(['deleting vertex ' num2str([v1])]);
			A(v1,:)=[];A(:,v1)=[];
			p=chrompoly(A);
			bp=[1 -1];tp=1;
			len=size(p,2);m3a=[];
			for i=1:size(p,2)
				m3a=[0 m3a] + p(len+1-i)*tp;
				tp=conv(tp,bp);
			end;
			out=[m3a 0];
		else
    
			if e < n*(n-1)/2,
				%     %choose 2 random vertices
				%     p=find(valseq>0);
				%     v1=p(ceil(rand(1)*size(p,2)));
				%     nb1=find(A(v1,:)==1);
				%     v2=nb1(ceil(rand(1)*valseq(v1)));

				nb1=find(A(v1,:)==1);
				nbv=valseq(nb1);
				v2=nb1(min(find(nbv==max(nbv))));

				%deletion
				%   disp(['deleting edge ' num2str([v1 v2])]);
				A1=A;
				A1(v1,v2)=0;A1(v2,v1)=0;
				m1=chrompoly(A1);
    
				%contraction
				%    disp(['contracting edge ' num2str([v1 v2])]);
				A2=A;
				A2(v1,:)=(A2(v2,:)+A2(v1,:))>0;
				A2(:,v1)=(A2(:,v2)+A2(:,v1))>0;
				A2(v1,v1)=0;A2(:,v2)=[];A2(v2,:)=[];
				m2=chrompoly(A2);
    
				s1=size(m1,2);
				s2=size(m2,2);
				pd=zeros(1,abs(s2-s1));
    
				if s1>s2,
					out=m1 - [pd m2];
				elseif s2>s1,
					out=[pd m1] - m2;
				else
					out= m1 - m2;
				end
			else
    
				nb1=find(A(v1,:)==0);
				nbv=valseq(nb1);
				v2=v1;
				while v2==v1, v2=nb1(min(find(nbv==max(nbv)))); end;

				%deletion
				%   disp(['adding edge ' num2str([v1 v2])]);
				A1=A;
				A1(v1,v2)=1;A1(v2,v1)=1;
				m1=chrompoly(A1);
    
				%contraction
				%    disp(['contracting edge ' num2str([v1 v2])]);
				A2=A;
				A2(v1,:)=(A2(v2,:)+A2(v1,:))>0;
				A2(:,v1)=(A2(:,v2)+A2(:,v1))>0;
				A2(v1,v1)=0;A2(:,v2)=[];A2(v2,:)=[];
				m2=chrompoly(A2);
    
				s1=size(m1,2);
				s2=size(m2,2);
				pd=zeros(1,abs(s2-s1));
    
				if s1>s2,
					out=m1 + [pd m2];
				elseif s2>s1,
					out=[pd m1] + m2;
				else
					out= m1 + m2;
				end

			end; %del-cont  / add-id
		end;
	end;
end; %if e==0



%%% OLD VERSION
% 
% n = nv(g);
% m = ne(g);
% 
% if m==0
%     p = zeros(1,n+1);
%     p(1) = 1;
%     return
% end
% 
% elist = edges(g);
% u = elist(1,1);
% v = elist(1,2);
% 
% h = graph;
% copy(h,g);
% 
% delete(h,u,v);
% p1 = chrompoly(h);
% contract(h,u,v);
% p2 = chrompoly(h);
% 
% p = p1 - [0,p2];
% 
% free(h);
