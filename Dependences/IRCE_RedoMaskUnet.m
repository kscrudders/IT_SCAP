function Base_label_ROIs = IRCE_RedoMaskUnet(Specific_ROIs, roi_corners, Ch1_corr_IRM, IRM_LUT, Save_individual_acq_dir, net)
% IRCE_RedoMaskUnet
% Manual mask cleanup per ROI using a key-driven editor (WindowKeyPressFcn).
% Keys: k=subtract, l=approve+next, ;=back, a=approve remaining & finish ROI
%
% Notes
% Base IRCE_MaskSeperation, with the initial label image replace with a
% unet segmentation output

    %---------------------------------------------------------%
    % Setup the mask for separating cell ROIs that are incorrectly connected
    %---------------------------------------------------------%
    subPath = fullfile(Save_individual_acq_dir, 'globalAnnotMaskSub.mat'); % path to subtracted mask annotations
    addPath = fullfile(Save_individual_acq_dir, 'globalAnnotMaskAdd.mat'); % path to additive mask annotations
    
    if exist(subPath,'file') && exist(addPath,'file')
        load(subPath,'globalAnnotMaskSub');
        load(addPath,'globalAnnotMaskAdd');
    else
        error('Must generate masks before editting them. Run Section_03b')
    end

    Base_label_ROIs = cell(size(roi_corners,1), 1);
    if isfile(fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat'))
        S2 = load(fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat'), 'Base_label_ROIs');
        if isfield(S2,'Base_label_ROIs')
            Base_label_ROIs = S2.Base_label_ROIs;
        end
    end


    % Loop over manually selected ROIs
    for i = Specific_ROIs
        fprintf('***** Current ROI = %03d *****\n', i);

        x = roi_corners{i,1}(:,1);
        y = roi_corners{i,1}(:,2);

        rr = min(y):max(y);                   % rows in ROI box
        cc = min(x):max(x);                   % cols in ROI box

        % Extract ROI stack and its auto segmentation
        imgROI = Ch1_corr_IRM(rr,cc,:);

        %---------------------------------------------------------%
        % Apply Unet
        %---------------------------------------------------------%
        Base_label = LF_apply_unet2stack(imgROI, net);

        %---------------------------------------------------------%
        % Clean up that basic Segmentation
        %---------------------------------------------------------%
        P = 202; % ~≥5 µm^2 @ 0.157 µm/px
        small_obj_removed = false(size(imgROI));
        for t = 1:size(imgROI,3)
            CC = bwconncomp(Base_label(:,:,t), 8);
            S  = regionprops(CC, 'Area');
            L  = labelmatrix(CC);
            keep = ismember(L, find([S.Area] >= P));
            small_obj_removed(:,:,t) = keep;
        end

        SE = strel('disk',6);
        Base_label = false(size(imgROI));
        for t = 1:size(imgROI,3)
            M = imdilate(small_obj_removed(:,:,t), SE);
            M = imfill(M, 'holes');
            M = imerode(M, SE);
            Base_label(:,:,t) = M;
        end

        % Apply global annotation maskes
        Base_label = (Base_label | globalAnnotMaskAdd(rr,cc,:)) .* globalAnnotMaskSub(rr,cc,:);

        %---------------------------------------------------------%
        % Manually Edit Masks
        %---------------------------------------------------------%
        [Base_label, globalAnnotMaskAdd(rr,cc,:), globalAnnotMaskSub(rr,cc,:)] = ...
            editMaskStack_viaKeys(imgROI, Base_label, globalAnnotMaskAdd(rr,cc,:), globalAnnotMaskSub(rr,cc,:), IRM_LUT);

        %---------------------------------------------------------%
        % Label the stack and save the edits
        %---------------------------------------------------------%
        Base_label_ROIs{i,1} = KLS_Label_fullStackTracking(Base_label);

        save(fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat'), 'Base_label_ROIs', '-v7.3');
        save(fullfile(Save_individual_acq_dir,'globalAnnotMaskSub.mat'), 'globalAnnotMaskSub', '-v7.3')
        save(fullfile(Save_individual_acq_dir,'globalAnnotMaskAdd.mat'), 'globalAnnotMaskAdd', '-v7.3')
    end

    clc;
    close all;

    %==================== Local functions (shared ws) ========================%
    function [Base_label, globalAnnotMaskAdd, globalAnnotMaskSub] = editMaskStack_viaKeys(imgROI, Base_label, globalAnnotMaskAdd, globalAnnotMaskSub, IRM_LUT)
    % Key-driven ROI mask editor with a main event loop that keeps the figure alive
    % until the user explicitly finishes the ROI.
    %
    % Keys: j=add, k=subtract, l=approve & next, ;=back one frame, a=approve remaining & finish
    
        [H,W,T] = size(imgROI);
    
        % ------------------- session state ------------------- %
        ii = 1;                          % current frame index
        initial_valid_seen = false;      % have we seen any valid interior boundary yet?
        approve_remaining = false;       % if true, end session
        finished = false;                % loop exit flag

        new_boundaries = (Base_label(:,:,1) | globalAnnotMaskAdd(:,:,1)) .* globalAnnotMaskSub(:,:,1);

    
        % --------------------- UI setup ---------------------- %
        f = figure('Name','ROI mask editor (j add, k sub, l next, ; back, a finish)', ...
                   'NumberTitle','off', ...
                   'Color','k', ...
                   'WindowKeyPressFcn',@onKey, ...
                   'Visible','off');
    
        ax = axes('Parent',f);
        ax.Visible = 'off';
        axis(ax,'image');
    
        imh = imshow(imgROI(:,:,1), IRM_LUT, 'Parent',ax, 'InitialMagnification','fit');
        hold(ax,'on');
        olRed  = plot(ax, nan, nan, 'r:', 'LineWidth', 2);  % edge-touching boundaries
        olBlue = plot(ax, nan, nan, 'b-', 'LineWidth', 2);  % interior boundaries
        hold(ax,'off');
    
        f.Visible = 'on';

        % Disable all pointer/scroll/btn tracking while idle (waiting for keys)
        toggleTracking(f, ax, false);
    
        % ------------------ MAIN EVENT LOOP ------------------ %
        while isgraphics(f) && ~finished
            % renderFrame (1)
            renderFrame(false); % render current frame (no back-jump)

            toggleTracking(f, ax, false);

            if ii > T
                continue
            end
            uiwait(f); % onKey will call uiresume(f)
        end
    
        if isgraphics(f)
            close(f);
        end
    
        % ================= nested helpers ==================== %
    
        function onKey(~,evt)
            key = evt.Key;
            if strcmp(key,'unknown') && isfield(evt,'Character') && ~isempty(evt.Character) && evt.Character==';'
                key = 'semicolon';
            end
    
            switch key
                case 'j'  % ADD (only when not in quick subtract-only mode)
                    toggleTracking(f, ax, true);
                    c = onCleanup(@() toggleTracking(f, ax, false));  % auto-disable on exit
                    doDraw(true); % Add to the mask
    
                case 'k'  % SUBTRACT
                    toggleTracking(f, ax, true);
                    c = onCleanup(@() toggleTracking(f, ax, false));      % auto-disable on exit
                    doDraw(false); % Subtract from the mask
    
                case 'l'  % approve current frame & advance by one
                    approveCurrentAndAdvance();
    
                    if ii > T
                        finished = true;
                    end
    
                case 'semicolon'  % go back one frame
                    goPrevFrame();
    
                case 'a'  % approve remaining & finish
                    Base_label(:,:,ii) = new_boundaries;
                    approve_remaining = true;
                    finished = true;
    
                otherwise
                    % ignore
            end
    
            if isgraphics(f)
                uiresume(f);
            end
        end
    
        function doDraw(isAdd)
            if isAdd
                roi = drawfreehand(ax, 'Color','g', 'LineWidth',5, ...
                                   'Closed',true, 'DrawingArea','auto', 'Multiclick',false);
            else
                roi = drawfreehand(ax, 'Color','r', 'LineWidth',5, ...
                                   'Closed',false, 'DrawingArea','auto', 'Multiclick',true);
            end

            if ~isgraphics(roi)
                return
            end
    
            SEline = strel('disk',3);
            M = createMask(roi, imh);
            M = imdilate(M, SEline);
            delete(roi);
    
            if isAdd
                globalAnnotMaskAdd(:,:,ii) = globalAnnotMaskAdd(:,:,ii) | M;
                globalAnnotMaskSub(:,:,ii) = globalAnnotMaskSub(:,:,ii) | globalAnnotMaskAdd(:,:,ii); % keep masks consistent

                new_boundaries = new_boundaries | M;
            else
                globalAnnotMaskSub(:,:,ii) = globalAnnotMaskSub(:,:,ii) & ~M;
                globalAnnotMaskAdd(:,:,ii) = globalAnnotMaskAdd(:,:,ii) & ~M; % keep masks consistent

                new_boundaries = new_boundaries & ~M;
            end
    
            updateOverlay();
            updateTitle();
        end
    
        function approveCurrentAndAdvance()
            Base_label(:,:,ii) = new_boundaries;
            ii = ii + 1;
        end
    
        function goPrevFrame()
            ii = max(1, ii - 1);
            isGoingBack = true;
        end
    
        function renderFrame(isGoingBack)
            % renderFrame (1) - keep track of ii as that is the current frame
            if nargin < 1
                isGoingBack = false;
            end
    
            if ii > T || approve_remaining
                finished = true;
                return
            end

            imh.CData = imgROI(:,:,ii);
            new_boundaries = (Base_label(:,:,ii) | globalAnnotMaskAdd(:,:,ii)) .* globalAnnotMaskSub(:,:,ii);


            if ii > T || approve_remaining
                finished = true;
                return
            end

            updateOverlay();
            updateTitle();
            drawnow limitrate;
        end
    
        function updateOverlay()
            [xr,yr,xb,yb] = prepBoundaryColor(new_boundaries, H, W);
            set(olRed,  'XData', xr, 'YData', yr);
            set(olBlue, 'XData', xb, 'YData', yb);
        end
    
        function updateTitle()
            modeTxt = "edit: add (j), subtract (k)";
            title(ax, sprintf("ROI %d | Frame %d/%d | %s | next (l), back (;), finish (a)", i, ii, T, modeTxt), ...
                  'Color',[0.9 0.9 0.9], 'FontWeight','bold');
        end
    end

    %================== Shared small helpers (local) ========================%
    function [xr,yr,xb,yb] = prepBoundaryColor(BW, rows, cols)
        % xr, yr - x and y values of masks that touch the edge, color 'red'
        % xb, yb - x and y values of masks that do not touch the edge,
            % color 'blue'
        B = bwboundaries(BW,'noholes');
        xr = [];
        yr = [];
        xb = [];
        yb = [];
        for k = 1:numel(B)
            boundary = B{k};
            touches = any(boundary(:,1)==1 | boundary(:,1)==rows | ...
                          boundary(:,2)==1 | boundary(:,2)==cols);
            if touches
                xr = [xr; boundary(:,2); NaN]; %#ok<AGROW>
                yr = [yr; boundary(:,1); NaN];
            else
                xb = [xb; boundary(:,2); NaN];
                yb = [yb; boundary(:,1); NaN];
            end
        end
    end
end

function toggleTracking(fig, ax, onoff)
    if onoff   % enable minimal interactions (for ROI drawing)
        enableDefaultInteractivity(ax);
        try
            iptPointerManager(fig,'enable'); 
        end
        fig.WindowButtonMotionFcn = [];  % let ROI manage motion
    else       % disable everything (idle waiting for keys)
        disableDefaultInteractivity(ax);
        datacursormode(fig,'off');
        fig.WindowButtonMotionFcn = [];
        fig.WindowScrollWheelFcn  = [];
        fig.WindowButtonDownFcn   = [];
        try
            iptPointerManager(fig,'disable'); 
        end
    end
end

function outputStack = LF_apply_unet2stack(inputStack, net)
    fillValue = 1; % Illumintation corrected IRM data has a bkgd value of 1

    % Get dimensions
    [x,y,t] = size(inputStack); % assume x,y,t order
    

    % Case 1: smaller than 512
    if x <= 512 && y <= 512
        padded = ones(512,512,t) * fillValue;
        padded(1:x,1:y,:,:) = inputStack;

        outputStack = false(size(padded));

        for i = 1:t
            pred = semanticseg(padded(:,:,i), net); 
            outputStack(:,:,i) = pred == 'cell';
        end

        outputStack = outputStack(1:x,1:y,:,:);

    % Case 2: larger than 512
    else
        xSplits = ceil(x/512);
        ySplits = ceil(y/512);
        xEdges = round(linspace(1,x+1,xSplits+1));
        yEdges = round(linspace(1,y+1,ySplits+1));

        outputStack = false(x,y,t);

        for ix = 1:xSplits
            for iy = 1:ySplits
                xs = xEdges(ix):xEdges(ix+1)-1;
                ys = yEdges(iy):yEdges(iy+1)-1;
                tile = inputStack(xs,ys,:);

                % pad each tile to 512×512
                tx = numel(xs);
                ty = numel(ys);
                padded = ones(512,512,t) * fillValue;
                padded(1:tx,1:ty,:) = tile;

                % run UNET slice by slice
                tMask = false(size(padded));
                for i = 1:t
                    pred = semanticseg(padded(:,:,i),net);
                    tMask(:,:,i) = pred == 'cell';
                end

                % crop back
                outputStack(xs,ys,:) = tMask(1:tx,1:ty,:);
            end
        end
    end
end