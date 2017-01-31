% Example script and functions for fitting a Gaussian to a set of measured data
%
% USAGE:
%	> gaussFitEx
%
%
% Christopher Smeenk
% January 30, 2017

function gaussFitEx()
% measured data: update these values
% ------------------------------------------------
positiona = [7.5,7.2,7.,6.8,6.6,6.4,6.3,6.2,6.0,5.8,5.6,5.4,5.2,5.]; % x
Ena = [3,8,15,26,38,47,50,50,47,39,25,15,7,3]; % y
Ena_err = ones(size(Ena)); % error in y
% ------------------------------------------------
% end of measured data


% centre the spots on zero for plotting
positiona = positiona - mean(positiona);

% plot the data
hn1 = errorbar(positiona, Ena, Ena_err, '-o');
xlabel('relative position (mm)');
ylabel('pulse energy (uJ)');
title('Wedge XF-532/266 spot size at 532 nm');
hold on;

% fit the Gaussian
% initial guess for parameters (centre,sigma,amplitude,offset)
iguess = [0,0.5,50,0];

% the next line calls the optimization function fminsearch to minimize the value of expfit by varying the parameter values in the vector p. p is a dummy variable.
fitp = fminsearch(@(p) expfit(p,positiona,Ena,Ena_err), iguess);

% plot the fit
xfit = linspace(positiona(1), positiona(end));
yfit = gaussFunc(xfit, fitp(1), fitp(2), fitp(3), fitp(4));
plot(xfit, yfit, 'r', 'LineWidth',2);
hold off;

% End of script
% ------------------------------------------------
end

% Functions to assist with fitting
% Gaussian function
function y = gaussFunc(x, x0, sigma, A, offset)
     %gaussFunc(x, x0, sigma, A):
     
     %OUTPUT - vector of values y evaluated at locations in x
	
     y = A*exp(-(x - x0).^2/sigma.^2) + offset;
end

% Chi-Squared function. This is the objective function to minimize by varying the gaussFunc parameters.
function chisqr = expfit(pars, xdata, ydata, errdata)
	%Fit a Gaussian using chi-squared minimization function.
	
	%FIT PARAMETERS:
	%pars(1) - x0
	%pars(2) - sigma
	%pars(3) - Amplitude
	%pars(4) - offset
	
	chisqr = sum( (ydata - gaussFunc(xdata, pars(1), pars(2), pars(3), pars(4)) ).^2/errdata.^2 );
end
