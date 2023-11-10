function [responses, psths] = combineResponsesAndPsths(cellArray)
% This function takes in a cell array that has the responses and psths from 
% different sessions and combines them (stacks them on top of each other)
% 
% INPUTS:
%     1) Cell Array 
%         - first column is responses (as produced by screeningScript))
%         - second column is psths (1x3 cell of raster, psth, times)
%             OR
%     2)Cell Array where first column is order and columns 2 and 3 are resp/psths respectively
%         
% OUTPUTS:
%     1) Responses (all cells)
%     2) Psths (all cells)
%     
% Note: It is important to confirm that you are only combining cells that saw the same stimuli
% vwadia Nov 2021

% Initial warning
disp('Did all these cells see the same stimuli?');
assert(size(cellArray, 2) > 1);

if size(cellArray, 2) > 2
    resp_col = 2; 
    psths_col = 3;
else
    resp_col = 1;
    psths_col = 2;
end
responses = {};
psths = {};
for row = 1:size(cellArray, 1)
    if size(cellArray{row, resp_col}, 2) > 3 % if basic method for resp lat computation
        responses = [responses; cellArray{row, resp_col}(:, 1:3)];
    else
        responses = [responses; cellArray{row, resp_col}];
    end
    psths = [psths; cellArray{row, psths_col}];

end