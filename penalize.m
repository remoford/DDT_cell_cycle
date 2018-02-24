function sum = penalize(objfun, numParams, params, parambounds)

outofbounds = 0;
for i=1:numParams
    if params(i) < parambounds(i,1) || params(i) > parambounds(i,2)
        outofbounds=1;
    end
end

if outofbounds == 1
    fprintf("out of bounds, penalizing!\n");
    sum = 0;
else
    sum = objfun(params);  
end


end