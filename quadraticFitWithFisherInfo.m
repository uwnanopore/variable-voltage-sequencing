function [params, fisher_info, total_entropy] = quadraticFitWithFisherInfo(x,y,varargin)
noise_uncertainty_multiplier = 1;
for ca = 1:length(varargin)
    switch upper(char(varargin{ca}))
        case 'MORENOISEUNCERTAINTY'
            noise_uncertainty_multiplier = 10;
            if length(varargin) > ca && isnumeric(varargin{ca+1})
                noise_uncertainty_multiplier = varargin{ca+1};
            end
    end
end
%perform the fit and find the noise calculated as the mean squared error
[x, y] = fitprep(x,y);
ft = fit(x,y,'poly2');
noise = sqrt(sum((ft(x) - y).^2)/length(x));
params = [ft.p1, ft.p2, ft.p3, noise];

%total relative entropy
total_entropy = length(x)*(log(2*pi*noise^2) + 1)/2;


%the fisher information is the hessian of the relative entropy
fisher_info = zeros(4);
fisher_info(1,1) = sum(x.^4)/noise^2;
fisher_info(1,2) = sum(x.^3)/noise^2;
fisher_info(2,1) = fisher_info(1,2);
fisher_info(2,2) = sum(x.^2)/noise^2;
fisher_info(3,1) = fisher_info(2,2);
fisher_info(1,3) = fisher_info(2,2);
fisher_info(2,3) = sum(x)/noise^2;
fisher_info(3,2) = fisher_info(2,3);
fisher_info(3,3) = length(x)/noise^2;
fisher_info(4,4) = 2*fisher_info(3,3);
fisher_info(4,1) = 2*sum((y - ft(x)).*x.^2)/noise^3;
fisher_info(1,4) = fisher_info(4,1);
fisher_info(4,2) = 2*sum((y - ft(x)).*x)/noise^3;
fisher_info(2,4) = fisher_info(4,2);
fisher_info(4,3) = 2*sum((y - ft(x)))/noise^3;
fisher_info(3,4) = fisher_info(4,3);

fisher_info(4,:) = fisher_info(4,:)/noise_uncertainty_multiplier;
fisher_info(:,4) = fisher_info(:,4)/noise_uncertainty_multiplier;

end