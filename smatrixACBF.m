function s = smatrixACBF(x1,K1,x2,K2,varargin)

pfbad = -inf;
selfalignment = false;
lookback = 5;

for ca = 1:length(varargin)
    switch upper(char(varargin{ca}))
        case 'PFBAD'
            pfbad = varargin{ca+1};
        case 'SELFALIGNMENT'
            lookback = varargin{ca+1};
            selfalignment = lookback > 0;
    end
end

if pfbad == -inf
    pfbad = -1e100;
end

% disp('computing score matrix')

%get the diagonals of the matrices for approximate bad feature removal
D1 = cell2mat(cellfun( @(x) diag(x), K1, 'uniformoutput', false));
D2 = cell2mat(cellfun( @(x) diag(x), K2, 'uniformoutput', false));


%if we are doing a self alignment
if selfalignment
    
    %preallocate the score matrix
    s = -inf(size(x1, 2), lookback+1);

    %loop over all levels
    for c1 = 1:size(x1, 2)
        
        %calculate matches for levels as far as the lookback
        for c2 = max(1, c1-lookback):c1
            
            %the difference in vectors squared
            dx_sq = (x1(:,c1) - x2(:,c2)).^2;
            
            %the combined variance
            V = (1./D1(:,c1) + 1./D2(:,c2));
            
            %compute the approximate score assuming no correlations,
            %presume this to be a good way to isolate "bad" features
            good_components = pfbad < -0.5*( log(2*pi*V) + dx_sq./V );
            
            %calculate the true score including the complete covariance
            %matrices, marginalizing the bad features.
            this_K1 = K1{c1}(good_components,good_components);
            this_K2 = K2{c2}(good_components,good_components);
            this_x1 = x1(good_components, c1);
            this_x2 = x2(good_components, c2);
            num_good = sum(good_components);
            Kx = this_K1*this_x1 + this_K2*this_x2;
            s(c1,c2-c1+lookback+1) = (size(x1, 1)-num_good)*pfbad + ...
                0.5*(  -num_good*log(2*pi) + sum(log(eig(this_K1))) + sum( log(eig(this_K2))) - sum(log(eig(this_K1+this_K2))) ...
                - this_x1'*this_K1*this_x1 - this_x2'*this_K2*this_x2 + (Kx'/(this_K1 + this_K2))*Kx);
        end
        
    end
    
    
else
    
    %preallocate the score matrix
    s = zeros(size(x1, 2), size(x2, 2));
    
    %loop over every combination of levels
    for c1 = 1:size(x1,2)
%         disp(c1)
%         disp([num2str(c1) '/1296'])
        
        for c2 = 1:size(x2,2)
            
            %the difference in vectors squared
            dx_sq = (x1(:,c1) - x2(:,c2)).^2;
            
            %the combined variance
            V = (1./D1(:,c1) + 1./D2(:,c2));
            
            %compute the approximate score assuming no correlations,
            %presume this to be a good way to isolate "bad" features
            good_components = pfbad < -0.5*( log(2*pi*V) + dx_sq./V );
            
            %calculate the true score including the complete covariance
            %matrices, marginalizing the bad features.
            this_K1 = K1{c1}(good_components,good_components);
            this_K2 = K2{c2}(good_components,good_components);
            this_x1 = x1(good_components, c1);
            this_x2 = x2(good_components, c2);
            num_good = sum(good_components);
            Kx = this_K1*this_x1 + this_K2*this_x2;
            s(c1,c2) = (size(x1, 1)-num_good)*pfbad + ...
                0.5*(  -num_good*log(2*pi) + sum(log(eig(this_K1))) + sum(log(eig(this_K2))) - sum(log(eig(this_K1+this_K2))) ...
                - this_x1'*this_K1*this_x1 - this_x2'*this_K2*this_x2 + (Kx'/(this_K1 + this_K2))*Kx);
            
            
        end
    end
    
end

end