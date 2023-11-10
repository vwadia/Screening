%% Compute the relevant features
%% Plotting the axes
% Need:
%     STA
%     Add synthetic_face_generator/ObjectSpace to the path (face/object screening)
%     Use stevens script to produce faces
%     make plot for that cell
%     run ScreeningScript until breakpoint
%   vwadia June2021
%%

% load params
if strcmp(taskStruct.subID, 'P73CS_Full')
    load([diskPath filesep 'ObjectSpace' filesep 'parameters_2k_synthetic_faces.mat'])
    params = params(1:667, :);
else 
    stimDir = '500Stimuli';
    layermat = 'fc6';
    load([diskPath filesep 'ObjectSpace' filesep stimDir filesep ['params_Alexnet_' layermat '_' stimDir '.mat']]);
end


%% Plot the preferred axis
% compute STA - there is a 1e5 scale difference between these two. 
% sta is bigger, linear_regression is smaller
% should not be an issue if one normalizes

% method = 'sta';
method = 'linear_regression';

for cellIndex = l(strctCells)
    if ~isempty(screeningData.responses{cellIndex, 2})
        [sta, ~] = analysis_STA(screeningData.responses{cellIndex, 1}, params, method);
        
        std_params = std(params); % std dev of each dimension
        if strcmp(method, 'linear_regression')
            
            if strcmp(taskStruct.subID, 'P73CS_Full')
                % shape-appearance dimensions are flipped between mine and Stevens code
                % why is this only for lin-reg?
                sta = fliplr(sta)./norm(sta); % normalize length
            else
                sta = sta./norm(sta);
            end
            all_sta_linReg(:, cellIndex) = sta;
            
        elseif strcmp(method, 'sta')
            sta = sta./norm(sta); % normalize length
            all_sta(:, cellIndex) = sta;
            
        end
        
        
        if strcmp(taskStruct.subID, 'P73CS_Full')          
            destPath = [basePath filesep 'rasters' filesep 'STA_Axis_Faces' filesep method];% filesep...
        elseif strcmp(taskStruct.subID, 'P73CS_ParamObj')
            destPath = [basePath filesep 'rasters' filesep 'STA_Axis_Objects' filesep method];% filesep...            
        end

        if ~exist(destPath)
            mkdir([destPath])
        end
            
        if strcmp(taskStruct.subID, 'P73CS_Full')
            % note that std_params is put in for the 'scale factor' -
            % because step n along the axis is n*stddev*STA, as the axis we
            % compute here of unit length 
            generate_face_50d_face_space(destPath,  'sta_with_range', sta, std_params, strctCells(cellIndex).Name);
        elseif strcmp(taskStruct.subID, 'P73CS_ParamObj')
            
            % Currently this projects the 19606 imagenet images into the
            % space built by the 1593 and then returns the nearest
            % neighbour of whatever point (feature vec) you give it)
            generate_objects_50d_obj_space(destPath, layer_resp, coeff1593, mu1593, im2all, 50, 'sta_with_range', sta, strctCells(cellIndex).Name);
        end
    end
end
 
%% plot the principal orthogonal axis




