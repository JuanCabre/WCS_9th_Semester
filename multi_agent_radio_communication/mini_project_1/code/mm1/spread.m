function out = spread(res,amp)

pdf = amp./sum(amp);
 
% The mean doppler shift
meanSpread = sum(res.*pdf); % Sum of each input in the pdf and res

varSum = sum(res.^2.*pdf);

varSpread = varSum - meanSpread^2;

out = sqrt(varSpread);

end
