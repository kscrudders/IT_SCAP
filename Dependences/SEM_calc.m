function sem = SEM_calc(vect, confLevel)
% SEM_calc - Standard Error of the Mean with specified confidence level
%
% Usage:
%   sem = SEM_calc(vect, confLevel)
%
% Inputs:
%   - vect: A vector or matrix. If a matrix, SEM is calculated for each column.
%           NaNs are ignored in both the mean and standard deviation calculations.
%
%   - confLevel [optional]: Confidence level (e.g., 0.95 for 95% confidence interval).
%                           Default is 0.95.
%
% Outputs:
%   - sem: Standard error of the mean multiplied by the critical value corresponding
%          to the specified confidence level.
%
% Notes:
%   - Uses the t-distribution to account for sample size.
%   - For large sample sizes, the t-distribution approximates the normal distribution.
%
% Example:
%   r = randn(1, 30);
%   S = SEM_calc(r, 0.99);
%   hist(r)
%   hold on
%   plot(mean(r), mean(ylim), 'r*')
%   plot([mean(r)-S, mean(r)+S], [mean(ylim), mean(ylim)], 'r-')
%   hold off
%

    % Input validation
    narginchk(1, 2);
    if nargin < 2 || isempty(confLevel)
        confLevel = 0.95; % Default confidence level
    end
    
    % Ensure 'vect' is a column vector or matrix
    if isvector(vect)
        vect = vect(:);
    end
    
    % Number of observations (excluding NaNs) in each column
    n = sum(~isnan(vect), 1);
    
    % Degrees of freedom
    df = n - 1;
    
    % Compute standard deviation (ignoring NaNs)
    stddev = std(vect,'omitnan');
    
    % Critical t-values for the specified confidence level
    tCritical = tinv(1 - confLevel / 2, df);
    
    % Handle cases where degrees of freedom are zero or negative
    tCritical(df <= 0) = NaN;
    
    % Compute standard error of the mean
    sem = (stddev ./ sqrt(n)) .* tCritical;
end