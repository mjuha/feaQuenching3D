function [fe,me,ke] = weakform(el,xe,de,ae,matNum)

global convectionLoad MAT fluxLoad

% dt = TS{1};
% alpha = TS{3};

% 1 point formula - degree of precision 1
gp =  [ 1/4, 1/4, 1/4;];
w = 1/6;
% gp = [0.5854101966249685  0.1381966011250105  0.1381966011250105; ...
%     0.1381966011250105  0.1381966011250105  0.1381966011250105; ...
%     0.1381966011250105  0.1381966011250105  0.5854101966249685; ...
%     0.1381966011250105  0.5854101966249685  0.1381966011250105];
% w = (1/6)*[0.25 0.25 0.25 0.25];

% get material properties
%prop = cell2mat(MAT(matNum));
prop = MAT(matNum);

ke = zeros(4,4);
me = zeros(4,4);
for i = 1:length(w)
    % stress-strain displacement matrix
    B = zeros(3,4);
    % loop over gauss points
    [N,dN,jac] = shape(gp(i,:),xe);
    % compute temperature
    T = N * de;
    % compute properties
    [ rho, k, cp ] = ComputeProperties( T, prop );
    for j=1:4 % loop over local nodes
        B(:,j) = dN(j,:);
    end
    ke = ke + B' * k * B * w(i) * jac;
    %
    me = me + N' * (rho*cp) * N * w(i) * jac;
end

% add contribution from temperature gradient
fe = ke * de;

if size(convectionLoad,1) > 0
    index = find(convectionLoad(:,1)==el,1); % 1 face
    flag = size(index,1);
    % compute side load
    if flag > 0
        [ke1,fe1] = computeSideLoad(index,xe,de,false);
        fe = fe + fe1;
        ke = ke + ke1;
    end
end

% now flux load
if size(fluxLoad,1) > 0
    index = find(fluxLoad(:,1)==el,1); % 1 face
    flag = size(index,1);
    % compute side load
    if flag > 0
        [~,fe1] = computeSideLoad(index,xe,de,true);
        fe = fe + fe1;
    end
end

%fe = fe - ke * de;
fe = fe + me*ae;
%
%me = me + alpha*dt*ke;

end
