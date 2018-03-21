function [ke,fe] = computeSideLoad(el,xe,de,isFluxLoad)

global convectionLoad fluxLoad HTableData isNBCTempDependent

% Gauss - Legendre rule
gp = [1/3, 1/3];
w = 0.5;
% gp = [  0.66666666666666666667  0.16666666666666666667; ...
%     0.16666666666666666667  0.66666666666666666667; ...
%     0.16666666666666666667  0.16666666666666666667];
% w = 0.5 * [0.33333333333333333333, ...
%     0.33333333333333333333, ...
%     0.33333333333333333333];


% compute residual: side loads
if ~isFluxLoad
    face = convectionLoad(el,2);
    if ~isNBCTempDependent
        h = convectionLoad(el,6); % coefficient
    end
    Ta = convectionLoad(el,7); % ambient temperature
else
    face = fluxLoad(el,2);
    q = fluxLoad(el,6); % heat flux
end


fe = zeros(4,1);
ke = zeros(4,4);
% loop over gauss points
for i = 1:length(w)
    % below we are going to use the same interchangably r,s and t.
    % But be aware that this is only for simplifying the program
    if face == 1 % local nodes 1-2-3
        r = gp(i,1);
        s = gp(i,2);
        % t = 0
        Nshape = [ 1-r-s, r, s, 0 ];
        N_r = [ -1, 1, 0, 0 ];
        N_s = [ -1, 0, 1, 0 ];
        % N_t = [0, 0, 0, 0];
    elseif face == 2 % local nodes 1-2-4
        r = gp(i,1);
        s = gp(i,2);
        % s = 0
        Nshape = [ 1-r-s, r, 0, s ];
        N_r = [ -1 ,1, 0, 0 ];
        % N_s = [0, 0, 0, 0];
        N_s = [ -1 ,0, 0, 1 ];
    elseif face == 3 % local nodes 1-3-4
        % r = 0
        r = gp(i,1);
        s = gp(i,2);
        Nshape = [ 1-r-s, 0, r, s ];
        % N_r = [0, 0, 0, 0];
        N_r = [ -1, 0, 1, 0 ];
        N_s = [ -1, 0, 0, 1 ];
    elseif face == 4 % local nodes 2-3-4
        % t = 1 - r - s
        r = gp(i,1);
        s = gp(i,2);
        Nshape = [ 0, r, s, 1-r-s ];
        N_r = [ 0, 1, 0, -1 ];
        N_s = [ 0, 0, 1, -1 ];
        % N_t = [ 0, 0, 0, 0 ];
    else
        error('Wrong face, check input!');
    end
    
    x_r = N_r * xe(:,1);
    x_s = N_s * xe(:,1);
    %
    y_r = N_r * xe(:,2);
    y_s = N_s * xe(:,2);
    %
    z_r = N_r * xe(:,3);
    z_s = N_s * xe(:,3);
    
    A = [x_r, y_r, z_r];
    B = [x_s, y_s, z_s];
    
    jac = norm(cross(A,B));
    
    % check jacobian
    if jac < 1.0e-14
        error('Negative jacobian, element too distorted!');
    end
    
    % surface temperature
    Ts = Nshape * de;
    
    if isNBCTempDependent
        h = interp1(HTableData(:,1), HTableData(:,2), Ts, 'linear', 0.0);
    end

    if ~isFluxLoad
        fe = fe - Nshape' * ( h * ( Ta - Ts ) ) * w(i) * jac;
        ke = ke + Nshape' * h * Nshape * w(i) * jac;
    else
        fe = fe - Nshape' * q * w(i) * jac;
        %ke = zeros(4,4);
    end
    
end
end