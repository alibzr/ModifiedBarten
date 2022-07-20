function output = sigmaEst(model,experimentalData)

summation = 0;

for i = 1:length(experimentalData)

summation = summation + (experimentalData(i) - model(i)).^2 ./ (length(experimentalData) - 2);

end

output = sqrt(summation);

end