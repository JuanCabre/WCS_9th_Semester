function out = spread(res,amp)

pdf = amp./sum(amp);


% calculate mean 
% The mean doppler shift
N = length(res);
spreadSum = sum(res.*pdf); % Sum of each input in the pdf and res
meanSpread = spreadSum; % Taking the mean of the 

varSum = sum(res.^2.*pdf);

varSpread = varSum - meanSpread^2;

out = varSpread;

end
