function [ rho, kappa, cp ] = ComputeProperties( T, prop )

% density
tmp = prop{1};
tempData = tmp{8};

Ta = tempData(1); % starting temperature for transformation of austenite
Tb = tempData(2); % starting temperature for transformation of perlite
    
if T > Ta % austenite
    rhoA = tmp{1};
    rho = interp1(rhoA(:,1), rhoA(:,2),T,'linear','extrap');
%     if isnan(rho)
%         error('Temperature is aout of range, please check input data!')
%     end
    cpA = tmp{3};
    cp = interp1(cpA(:,1), cpA(:,2),T,'linear','extrap');
%     if isnan(cp)
%         error('Temperature is aout of range, please check input data!')
%     end
    kA = tmp{5};
    kappa = interp1(kA(:,1), kA(:,2),T,'linear','extrap');
%     if isnan(kappa)
%         error('Temperature is aout of range, please check input data!')
%     end
elseif T > Tb
    rhoA = tmp{2};
    rho = interp1(rhoA(:,1), rhoA(:,2),T,'linear','extrap');
%     if isnan(rho)
%         error('Temperature is aout of range, please check input data!')
%     end
    cpA = tmp{4};
    cp = interp1(cpA(:,1), cpA(:,2),T,'linear','extrap');
%     if isnan(cp)
%         error('Temperature is aout of range, please check input data!')
%     end
    kA = tmp{6};
    kappa = interp1(kA(:,1), kA(:,2),T,'linear','extrap');
%     if isnan(kappa)
%         error('Temperature is aout of range, please check input data!')
%     end
else % martensite
    rhoA = tmp{1};
    rho = interp1(rhoA(:,1), rhoA(:,2),T,'linear','extrap');
%     if isnan(rho)
%         error('Temperature is aout of range, please check input data!')
%     end
    cpA = tmp{3};
    cp = interp1(cpA(:,1), cpA(:,2),T,'linear','extrap');
%     if isnan(cp)
%         error('Temperature is aout of range, please check input data!')
%     end
    kA = tmp{7};
    kappa = interp1(kA(:,1), kA(:,2),T,'linear','extrap');
%     if isnan(kappa)
%         error('Temperature is aout of range, please check input data!')
%     end
end


end

