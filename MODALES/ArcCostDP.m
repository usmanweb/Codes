% Source: D. P. Bertsekas, Dynamic Programming

function [route_st]=ArcCostDP(L,s,t,steps)

%% Inputs
% L = Cost Matrix
% s = min speed value 
% t = max speed value 
% steps = # of arcs

n = size(L,2);

L(find(L==0))=Inf;  % make all zero distances equal to infinity

for i=1:n
  J(steps,i) = L(i,t); 
  route(steps,i).path = [t];
end

% find min for every i: Jk(i)=min_j(L(i,j)+Jk+1(j))
for p=1:steps-1
  k=steps-p; % recurse backwards
  
  for i=1:n
    %fprintf('stage %2i, node %2i \n',k,i)
    [J(k,i),ind_j] =  min(L(i,:)+J(k+1,:));
    route(k,i).path = [ind_j, route(k+1,ind_j).path];
  end
  
end

[J_st,step_ind] = min(J(:,s));
route_st = [s, route(step_ind,s).path];
J=J(sort(1:n,'descend'),:);
route=route(sort(1:n,'descend'),:);