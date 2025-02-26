function [fitFunction, gof, fit_str, lifetime_tau] = KLS_Exponentialfit_and_plot(imageStack, plot_flag)
    % Corrects a 3D image stack for photobleaching using an exponential fit and plots the fit.
    %
    % Parameters:
    %   imageStack - 3D array of image slices (height x width x numSlices)
    %
    % Returns:
    %   fitFunction - function handle for the exponential fit
    
    if nargin < 2
        plot_flag = 0;
    end
    
    % Get the size of the image stack
    [~, ~, numSlices] = size(imageStack);

    % Calculate the mean intensity of each slice
    %meanIntensities = zeros(numSlices, 1);
    %for i = 1:numSlices
    %    meanIntensities(i) = mean(imageStack(:,:,i), 'all');
    %end
    
    y = squeeze(mean(imageStack,[1 2])); % Pull the mean value for the 

    % Fit the mean intensities to an exponential decay model with a constant term
    x = (1:numSlices)'-1;
    
    ft = fittype('a*exp(-b*x) + c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMaxChange = 0.01;
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.MaxIter = 1000;
    
    IG_A = max(y); % Initial guess for A
    
    y_at_taulifetime = max(y)-(max(y)*log(2));
    [~, index] = min(abs(y-y_at_taulifetime)); % find a very rough approximation for time at tau half-life
    IG_B = 1/index; % Initial guess for B
    
    IG_C = min(y); % Initial guess for C
    opts.StartPoint = [IG_A IG_B IG_C];
    opts.Robust = 'Bisquare';    
    [fitParams, gof] = fit(x, y, ft, opts);

    % Exponential fit function with constant term
    fitFunction = @(x) fitParams.a * exp(-fitParams.b * x) + fitParams.c;

    % Get the confidence intervals for the fit
    ci = confint(fitParams, 0.95);
    fitFunctionLower = @(x) ci(1,1) * exp(-ci(1,2) * x) + ci(1,3);
    fitFunctionUpper = @(x) ci(2,1) * exp(-ci(2,2) * x) + ci(2,3);
    
    if plot_flag == 1
        % Plot the mean intensities and the exponential fit
        plot(x, y, 'bo', 'DisplayName', 'Data');
        hold on;
        plot(x, fitFunction(x), 'r-', 'DisplayName', ...
            sprintf('Fit: y = %.0f * e^{(-%.4f * x)} + %.0f\nR^2 = %.3f', fitParams.a, fitParams.b, fitParams.c, gof.rsquare));
        
        % Shade the area between the confidence intervals
        xVals = [x; flipud(x)];
        yVals = [fitFunctionLower(x); flipud(fitFunctionUpper(x))];
        fill(xVals, yVals, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        plot(x, fitFunctionLower(x), 'r--', 'DisplayName', '95% CI Lower');
        plot(x, fitFunctionUpper(x), 'r--', 'HandleVisibility', 'off');

        xlabel('Frame');
        ylabel('Mean Intensity (AU)');
        title('Photobleaching Correction: Exponential Fit');
        legendHandle = legend('show');
        set(legendHandle, 'FontSize', 9, 'location', 'northeast'); % Set the legend font size to 9
        hold off;        
    end
    
    fit_str = [num2str(fitParams.a,3) '*exp^{(-' num2str(fitParams.b,3) '*x)} +' num2str(fitParams.c,3)];
    lifetime_tau = 1/fitParams.b;
end
