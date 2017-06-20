function [c] = OMP( y, D, s)
%OMP Realization of Orthogonal Matching Pursuit
% Solving sparse optimization with the form ||y-Dx||^2

dict_column=zeros(1,s);
res=y;

for i=1:s
    proj = abs(res*D)./sqrt(sum(D.^2));
    % pick the max influence projection column
    max_v=-1;
    for j=1:length(proj)
        if isempty(find(dict_column==j,1))==0
            continue;
        end
        if max_v< proj(j)
            max_v=proj(j);
            index= j;
        end
    end
    % record its index
    dict_column(i)=index;
    % solve LSE problem by Moore inv of the matrix
    x= y/(D(:,dict_column(1:i))');
    
    % calculate the residual
    res= y- (D(:,dict_column(1:i))*x')';
    %norm(res)
    if norm(res) < 1e-3
        % break in advance
        c=zeros(1,size(D,2));
        c(dict_column(1:i))=x';
        return 
    end
end

% extend the basis (add zeros)
c=zeros(1,size(D,2));
c(dict_column)=x';
end

