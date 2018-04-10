function feaQuenching3D(filename,varargin)
% feaHeat3D.m
% This is the main file that implement a tetrahedra
% element for solving transient heat transfer problems
%
% Author: Dr. Mario J. Juha
% Date: 11/09/2017
% Mechanical Engineering
% Universidad de La Sabana
% Chia -  Colombia
%
% Clear variables from workspace
%clearvars

global nel neq nzmax coordinates U elements nn LM irow icol ID TS isTimeDBC
global NLOPT Phase

% Specify file name
%filename = '\\Client\C$\Users\marioju\Documents\Work\valvula\example.inp';

% check if user have specified rregression tests
if nargin > 1
    error('Only accept filename!')
%     % extract the regression values
%     normDisp = varargin{1};
%     normRot = varargin{2};
%     regressionTol = varargin{3};
end

% read data
fprintf('************************\n')
fprintf('Reading input data\n')
fprintf('************************\n\n')
fflush(stdout);
outfile = readData(filename);

% at time = 0
Tolerance = NLOPT(3);
maxNumIter = NLOPT(4);
fprintf('\nComputing time %f\n\n',0)
fflush(stdout);
if isTimeDBC
    DBC_InTime(0)
end
dt = TS{1}; % delta time
% Newton start here
newtonIter = 1;
while ( newtonIter <= maxNumIter )
    fprintf('  Computing iteration %d\n',newtonIter)
    fflush(stdout);
    % ===========================
    % assembling stiffness matrix
    % ===========================
    K = zeros(1,nzmax);
    F = zeros(neq,1);
    % set counter to zero
    count = 0;
    for i=1:nel
        xe = coordinates(elements(i,2:5),:);
        de = U(1,elements(i,2:5));
        ae = U(2,elements(i,2:5));
        % get phases on element�s nodes
        pe = Phase(:,elements(i,2:5));
        matNum = elements(i,1); % element material number
        [fe,me,~] = weakform(i,xe,de',ae',pe,matNum);
        for k=1:4
            i_index = LM(k,i);
            if (i_index > 0)
                F(i_index) = F(i_index) - fe(k);
                for m=1:4
                    j_index = LM(m,i);
                    if (j_index > 0)
                        count = count + 1;
                        K(count) = me(k,m);
                    end
                end
            end
        end
    end
    fprintf('\n')
    fprintf('Solving system of equations\n')
    fprintf('\n')
    fflush(stdout);
    M = sparse(irow,icol,K,neq,neq);
    F = M\F;
    % assign solution
    for r=1:nn
        i_index = ID(r);
        if (i_index > 0)
            U(2,r) = U(2,r) + F(i_index);
        end
    end
    % compute residual norm using new values
    normF = computeF();
    %         fprintf('  Norm of residual old: %g\n', oldNormF)
    fprintf('  Norm of residual: %g\n\n', normF)
    fflush(stdout);
    if normF < Tolerance
        fprintf('  Solution converged!\n')
        fflush(stdout);
        break
    end
    newtonIter = newtonIter + 1;
end % newtonIter
if newtonIter > maxNumIter
    error('Does not converge after max number of iterations!')
end
% now compute phases
fprintf('\n')
fprintf('Computing phases\n\n')
fflush(stdout);
phaseTransformation(dt)
WriteVTKFile(outfile,0)

% fprintf('************************\n')
% fprintf('Assembling stiffness matrix and force vector\n')
% fprintf('************************\n\n')
fprintf('************************\n')
fprintf('Advancing in time\n')
fprintf('************************\n\n')
fflush(stdout);
%
t = 0;
tf = TS{2}; % final time
alpha = TS{3}; % method
nts = TS{4}; % ouput every nts
counter = 0;
count1 = 0;
%
% newtonOpt = NLOPT(1);
% normOpt = NLOPT(2);
% Tolerance = NLOPT(3);
% maxNumIter = NLOPT(4);
while ( t < tf )
    counter = counter + 1;
    t = t + dt;
    fprintf('\nComputing time %f\n\n',t)
    fflush(stdout);
    if isTimeDBC
        DBC_InTime(t)
    end
%     oldNormF = computeF();
    % predictor value
    dp = U(1,:) + (1-alpha)*dt*U(2,:);
    U(1,:) = dp;
    U(2,:) = 0;
    % Newton start here
    newtonIter = 1;
    while ( newtonIter <= maxNumIter )
        fprintf('  Computing iteration %d\n',newtonIter)
        fflush(stdout);
        % ===========================
        % assembling stiffness matrix
        % ===========================
        K = zeros(1,nzmax);
        F = zeros(neq,1);
        % set counter to zero
        count = 0;
        for i=1:nel
            xe = coordinates(elements(i,2:5),:);
            de = U(1,elements(i,2:5));
            ae = U(2,elements(i,2:5));
            % get phases on element�s nodes
            pe = Phase(:,elements(i,2:5));
            matNum = elements(i,1); % element material number
            [fe,me,ke] = weakform(i,xe,de',ae',pe,matNum);
            for k=1:4
                i_index = LM(k,i);
                if (i_index > 0)
                    F(i_index) = F(i_index) - fe(k);
                    for m=1:4
                        j_index = LM(m,i);
                        if (j_index > 0)
                            count = count + 1;
                            K(count) = me(k,m)+ alpha*dt*ke(k,m);
                        end
                    end
                end
            end
        end
        fprintf('\n')
        fprintf('Solving system of equations\n')
        fprintf('\n')
        fflush(stdout);
        M = sparse(irow,icol,K,neq,neq);
        F = M\F;
        % assign solution
        for r=1:nn
            i_index = ID(r);
            if (i_index > 0)
                U(2,r) = U(2,r) + F(i_index);
                % corrector
                U(1,r) = dp(r) + alpha*dt*U(2,r);
            end
        end
        newtonIter = newtonIter + 1;
        % compute residual norm using new values
        normF = computeF();
        fprintf('  Norm of residual: %g\n\n', normF)
        fflush(stdout);
        if normF < Tolerance
            fprintf('  Solution converged!\n')
            fflush(stdout);
            break
        end
    end % newtonIter
    % now compute phases
    fprintf('\n')
    fprintf('Computing phases\n\n')
    fflush(stdout);
    phaseTransformation(dt)
    if mod(counter,nts) == 0
        count1 = count1 + 1;
        WriteVTKFile(outfile,count1)
    end
    if newtonIter > maxNumIter
        error('Does not converge after max number of iterations!')
    end
end
% fprintf('************************\n')
% fprintf('Computing stresses and printing final mesh and results\n')
% fprintf('************************\n\n')
% 
% computeStressStrain
% 
% WriteVTKFile(outfile,1)
end
