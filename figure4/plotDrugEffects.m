function [x, y, A, B, C] = plotDrugEffects(data)
    % Function to convert 3-component vectors to (x, y) positions in a triangle
    % and plot the points.
    %
    % Input:
    %   data - Nx3 matrix where each row is [ gen_time, max_load,lag],
    %   already normalized 
    
    % Normalize each row so the scores sum to 1
    normalizedData = data ./ sum(data, 2);

    % Triangle vertices
    A = [0, 0];                      % Bottom-left
    B = [1, 0];                      % Bottom-right
    C = [0.5, sqrt(3)/2];            % Top vertex
    
    % Compute (x, y) positions
    x = normalizedData(:, 2) * B(1) + normalizedData(:, 3) * C(1);
    y = normalizedData(:, 3) * C(2); % Only the w component contributes to y
end
