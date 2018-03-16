function [ normF ] = computeF
global neq nel coordinates elements LM U

F = zeros(neq,1);
for i=1:nel
    xe = coordinates(elements(i,2:5),:);
    de = U(1,elements(i,2:5));
    ae = U(2,elements(i,2:5));
    matNum = elements(i,1); % element material number
    [fe,~,~] = weakform(i,xe,de',ae',matNum);
    for k=1:4
        i_index = LM(k,i);
        if (i_index > 0)
            F(i_index) = F(i_index) + fe(k);
        end
    end
end

normF = norm(F);

end

