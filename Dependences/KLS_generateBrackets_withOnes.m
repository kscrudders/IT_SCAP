function KLS_generateBrackets_withOnes(n)
    % Function to print to the command line the text with a dynamic number of brackets
    % Input: n - the number of brackets to generate
    % The indentation and formatting follow the pattern provided
    
    fprintf('labelsForTheTargetCell_Post = {\n');
    
    totalLines = 0;
    cycleLength = 10; % The pattern cycles every 12 lines
    for i = 1:n
        % Adjust indentation based on position in the cycle
        positionInCycle = mod(i-1, cycleLength) + 1;
        if positionInCycle <= 5
            indentLevel = 1;
        else
            indentLevel = 2;
        end

        indent = repmat('    ', 1, indentLevel);


        totalLines = totalLines + 1;

        % Check if we need to print the line count comment
        if mod(totalLines, 10) == 0
            % Print the line count comment
            fprintf('%s[1]; ... %%0%02d \n', indent, totalLines);
        else
            % Print the line
            fprintf('%s[1]; ...\n', indent);            
        end
    end

    fprintf('};\n');
end