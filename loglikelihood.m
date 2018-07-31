function sum = loglikelihood(data, pdf, numParams, params )

numData = length(data);

probabilities = [];
sum=0;
switch numParams
    case 0
        probabilities = pdf(data);
    case 1
        probabilities = pdf(data,params(1));
    case 2
        probabilities = pdf(data,params(1), params(2));
    case 3
        probabilities = pdf(data,params(1), params(2), params(3));
    case 4
        probabilities = pdf(data,params(1), params(2), params(3), params(4));
    case 5
        probabilities = pdf(data,params(1), params(2), params(3), params(4), params(5));
    case 6
        probabilities = pdf(data,params(1), params(2), params(3), params(4), params(5), params(6));
    otherwise
        %fprintf("ERROR: unspported number of arguments\n");
end
for i=1:numData
    sum = sum + log(probabilities(i));
end

end