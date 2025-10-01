function [roi_corners, previous_ROI_count] = IRCE_CellROIs(Ch1_corr_IRM, Save_individual_acq_dir, Add_more_ROIs_flag, IRM_thres, IRM_LUT, ROIs_dir)
% Keys:
%   j : begin drawing a rectangle (drag, release to finish)
%   k : commit current ROI and start another
%   l : commit current ROI (if any) and finish  <-- loop exits on 'l'
%   ; : undo last committed ROI
%   a : finish (commit current ROI if present)  <-- also exits loop
% Updated 20240822: KLS - Built in check for x/y positions outside image dimensions
% KLS update (keypress loop) — 2025-09-27

% ---------- projections ----------
imgMin = min(Ch1_corr_IRM,[],3);
if size(Ch1_corr_IRM,3) > 15
    imgMedMin = min(movmedian(Ch1_corr_IRM, 15, 3), [], 3);
else
    win = max(3, ceil(size(Ch1_corr_IRM,3)*0.1));
    imgMedMin = min(movmedian(Ch1_corr_IRM, win, 3), [], 3);
end
median_mask = imgMedMin < IRM_thres;


% ---------- load any existing ROIs ----------
roiFile = fullfile(Save_individual_acq_dir,'roi_corners.mat');
roi_corners = {};
if isfile(roiFile)
    S = load(roiFile,'roi_corners');
    if isfield(S,'roi_corners')
        roi_corners = S.roi_corners;
    end
end
previous_ROI_count = size(roi_corners,1);

% ---------- if no new ROIs requested AND file exists, just render final & export ----------
if ~Add_more_ROIs_flag && ~isempty(roi_corners)
    finalizeAndExport();
    return
end

% ---------- overlay boundaries (edge-touching = red, interior = blue) ----------
[rows, cols] = size(median_mask);
[B, ~] = bwboundaries(median_mask > 0, 'noholes');
[xr,yr,xb,yb] = collectBoundaryPolylines(B, rows, cols);

% ---------- initial display ----------
f = figure('Name','Select cell ROIs (k - add ROI, l - finish, ; - undo)', ...
           'NumberTitle','off', ...
           'Color','k', ...
           'WindowKeyPressFcn',@onKey);
g = gcf;
g.WindowState = 'maximized';
ax = axes('Parent',f);
axis(ax,'image');
ax.Visible = 'off';
imshow(imgMin, IRM_LUT, 'Parent',ax, 'InitialMagnification','fit');
hold(ax,'on');
plot(ax, xr, yr, 'r:', 'LineWidth', 1);
plot(ax, xb, yb, 'b:', 'LineWidth', 1);
hold(ax,'off');
drawnow;



% ---------- render any existing rectangles (yellow) ----------
rectH = gobjects(0);
labelH = gobjects(0);
for i = 1:previous_ROI_count
    [x,y] = cornersToXY(roi_corners{i,1});
    [x,y] = stayWithinImage(size(imgMin), x, y);
    [pos, cx, cy] = cornersToRect([x,y]);
    rectH(i,1)  = drawrectangle(ax,'Position',pos,'Color','y','LineWidth',2,'FaceAlpha',0.03,'InteractionsAllowed','none');
    labelH(i,1) = text(ax, cx, cy, sprintf('%03d',i), 'Color','c','FontSize',14, ...
                       'HorizontalAlignment','center','FontName','Trebuchet MS','FontWeight','bold');
end



% ---------- interactive session state ----------
curRect   = [];                          % active drawrectangle handle
activeStarted = false;                   % actively drawing?
nextIndex = size(roi_corners,1) + 1;    % next ROI number
lastKey   = '';                          % last key pressed (normalized)
editTargetIdx = 0;     % 0 = append mode; >0 means we are editing that ROI index


% ---------- MAIN EVENT LOOP ---------- %
while isgraphics(f) && ~strcmp(lastKey,'l')
    figure(f)
    updateTitle();
    toggleTracking(f,ax,false);
    uiwait(f);                           % onKey calls uiresume to return here
    % loop continues until lastKey becomes 'l' (or fig closed)
end

% Close the interactor fig if still open
if isgraphics(f)
    close(f);
end

% ---------- show final numbers on the min projection and export ----------
finalizeAndExport()

% ======== nested functions (UI + helpers) ========

function onKey(~,evt)
    % Normalize key; MATLAB sends char like 'j','k','l','semicolon'
    key = evt.Key;
    if strcmp(key,'unknown')
        if isfield(evt,'Character') && ~isempty(evt.Character) && evt.Character==';'
            key = 'semicolon';
        end
    end

    switch key
        case 'k'      % commit & start another
            if commitCurrent()
                toggleTracking(f,ax,true);
                startDrawing();
            elseif isempty(curRect) || ~isgraphics(curRect)
                toggleTracking(f,ax,true);
                startDrawing();
            end

        case 'semicolon'   % undo last committed
            toggleTracking(f,ax,true);
            undoLast();
            toggleTracking(f,ax,false);

        case {'quote'}  % apostrophe pressed
            promptModifyOrRemove();

        case 'l'      % commit and finish
            toggleTracking(f,ax,true);
            commitCurrent();
            toggleTracking(f,ax,false);
            lastKey = 'l';

        otherwise
            % ignore other keys; keep looping
    end

    % Always return control to the main loop after any key
    if isgraphics(f)
        uiresume(f);
    end
end

function startDrawing()
    if isempty(curRect) || ~isgraphics(curRect)
        curRect = drawrectangle(ax,'Color','y','LineWidth',3,'FaceAlpha',0.05);
        activeStarted = true;
    end
    updateTitle();
end

function ok = commitCurrent()
    ok = false;
    if ~isempty(curRect) && isgraphics(curRect)
        pos = curRect.Position;            % [x y w h]
        if any(pos(3:4) <= 0)
            delete(curRect);
            curRect = [];
            activeStarted = false;
            updateTitle();
            return
        end

        % Clamp to image extent and round to pixel corners
        [pos, xy] = clampRectToImage(pos, size(imgMin));

        if editTargetIdx > 0
            % --------- REPLACE existing ROI ---------
            target = editTargetIdx;

            % Remove old persistent graphics for that ROI
            if target <= numel(rectH) && isvalidHandle(rectH(target))
                delete(rectH(target));
                rectH(target) = gobjects(1);
            end
            if target <= numel(labelH) && isvalidHandle(labelH(target))
                delete(labelH(target));
                labelH(target) = gobjects(1);
            end

            % Update stored corners in-place
            roi_corners{target,1} = xy;

            % Redraw persistent rectangle + label with same index
            [cx,cy] = rectCenter(pos);
            rectH(target,1)  = drawrectangle(ax,'Position',pos,'Color','y','LineWidth',2,'FaceAlpha',0.03,'InteractionsAllowed','none');
            labelH(target,1) = text(ax, cx, cy, sprintf('%03d',target), 'Color','c','FontSize',14, ...
                                     'HorizontalAlignment','center','FontName','Trebuchet MS','FontWeight','bold');

            % Clear edit mode
            editTargetIdx = 0;

        else
            % --------- APPEND new ROI (normal path) ---------
            roi_corners{nextIndex,1} = xy;     %#ok<AGROW>

            [cx,cy] = rectCenter(pos);
            rectH(nextIndex,1)  = drawrectangle(ax,'Position',pos,'Color','y','LineWidth',2,'FaceAlpha',0.03,'InteractionsAllowed','none');
            labelH(nextIndex,1) = text(ax, cx, cy, sprintf('%03d',nextIndex), 'Color','c','FontSize',14, ...
                                        'HorizontalAlignment','center','FontName','Trebuchet MS','FontWeight','bold');
            nextIndex = nextIndex + 1;
        end

        % Dispose editor rect
        delete(curRect);
        curRect = [];
        activeStarted = false;
        ok = true;
        updateTitle();

        % Persist safely
        try
            save(roiFile,'roi_corners','-v7.3');
        catch
            warning('Could not save ROI file:');
        end

        % Back to idle tracking
        try
            toggleTracking(f,ax,false);
        catch
        end
    end
end

function undoLast()
    % Cancel active editor first
    if ~isempty(curRect) && isgraphics(curRect)
        delete(curRect);
        curRect = [];
        activeStarted = false;
        updateTitle();
        return
    end

    if nextIndex > 1
        last = nextIndex - 1;

        if last >= 1 && last <= numel(rectH) && isvalidHandle(rectH(last))
            delete(rectH(last));
            rectH(last) = gobjects(1);
        end
        if last >= 1 && last <= numel(labelH) && isvalidHandle(labelH(last))
            delete(labelH(last));
            labelH(last) = gobjects(1);
        end
        if ~isempty(roi_corners)
            roi_corners(last,:) = [];
        end

        nextIndex = last;

        try
            save(roiFile,'roi_corners','-v7.3');
        catch
            warning('Could not save ROI file:');
        end
    end

    if editTargetIdx > 0
        % Cancel in-progress edit
        if ~isempty(curRect) && isgraphics(curRect)
            delete(curRect);
        end
        curRect = [];
        activeStarted = false;
        editTargetIdx = 0;
        updateTitle();
        return
    end

    updateTitle();
end

function updateTitle()
    mode = 'idle';
    if activeStarted
        mode = 'drawing (drag to size; k commit, l finish)';
    end
    title(ax, sprintf("ROIs: %d | k - add ROI, l - finish, ; - undo, ' - modify ROI — %s", ...
          max(0,nextIndex-1), mode), 'Color',[0.9 0.9 0.9],'FontWeight','bold');
    drawnow limitrate;
end

function promptModifyOrRemove()
    if size(roi_corners,1) == 0
        warndlg('No ROIs to modify/remove yet.','ROI Edit','modal');
        return
    end

    % Ask which ROI number to act on
    defNum = num2str(max(1, min(size(roi_corners,1), nextIndex-1)));
    ansNum = inputdlg( ...
        sprintf('Enter ROI # to modify (1–%d):', size(roi_corners,1)), ...
        'Select ROI', 1, {defNum});
    if isempty(ansNum)
        return
    end
    idx = str2double(ansNum{1});
    if ~(isfinite(idx) && idx>=1 && idx<=size(roi_corners,1) && idx==floor(idx))
        warndlg('Invalid ROI number.','ROI Edit','modal');
        return
    end
    idx = double(idx);

    %{
    % Ask action: remove (n) or modify (m)
    ansAct = inputdlg( ...
        sprintf('Action for ROI %d? Enter "n" = remove, "m" = modify:', idx), ...
        'ROI Action', 1, {'m'});
    if isempty(ansAct)
        return
    end
    %}
    % Ask action with custom button labels
    choice = questdlg( ...
        sprintf('ROI %d: What would you like to do?', idx), ...
        'ROI Action', ...
        'Remove ROI', 'Modify ROI', 'Modify ROI');   % default = Modify
    
    % Map to act: 'n' = remove, 'm' = modify
    if isempty(choice)
        % user cancelled — keep focus and bail
        try
            uicontrol(focusCtl); 
        catch
            figure(f); 
        end
        return
    end
    
    switch choice
        case 'Remove ROI'
            act = 'n';
        case 'Modify ROI'
            act = 'm';
        otherwise
            % safety fallback
            act = 'm';
    end


    switch act
        case 'n'  % remove (with filesystem cleanup + renumber)
        
            % --- If there are any files/folders, confirm destructive action ---
            hasAny = ~isempty(dir(fullfile(ROIs_dir,'*')));
            if hasAny
                msg = sprintf(['Delete ROI %03d?\n\n' ...
                               'This will permanently delete folder "%s" (if it exists) and its contents,\n' ...
                               'then shift higher-numbered folders down by one (Cell_003→Cell_002, etc.).\n\n' ...
                               'It will also remove entry %d from Base_label_ROIs.mat (if present).'], ...
                               idx, sprintf('Cell_%03d',idx), idx);
                choice = questdlg(msg, 'Confirm ROI delete', 'Delete', 'Cancel', 'Delete');
                if ~strcmp(choice,'Delete')
                    % Keep focus and bail without changing anything
                    try 
                        uicontrol(focusCtl) 
                    catch
                        figure(f); 
                    end
                    return
                end
            end
        
            % --- Delete Cell_### and renumber higher indexes on disk ---
            try
                deleteAndRenumberCellFolders(idx, ROIs_dir);
            catch
                warning('Failed to delete/renumber ROI folders:');
            end
        
            % --- Remove the same index from Base_label_ROIs.mat (if present) ---
            try
                blrFile = fullfile(Save_individual_acq_dir,'Base_label_ROIs.mat');
                if isfile(blrFile)
                    Sblr = load(blrFile,'Base_label_ROIs');
                    if isfield(Sblr,'Base_label_ROIs') && iscell(Sblr.Base_label_ROIs)
                        if idx>=1 && idx<=numel(Sblr.Base_label_ROIs)
                            Sblr.Base_label_ROIs(idx,:) = [];
                            save(blrFile,'-struct','Sblr','-v7.3');
                        end
                    end
                end
            catch
                warning('Failed to update Base_label_ROIs.mat:');
            end
        
            % --- Remove from current session (graphics + roi_corners) ---
            removeROI(idx);
        
            % Persist roi_corners after removal
            try
                save(roiFile,'roi_corners','-v7.3');
            catch
                warning('Could not save ROI file:');
            end
        
            % Return focus for more keys
            try
                uicontrol(focusCtl);
            catch
                figure(f);
            end
        %{
        case 'n'  % remove
            removeROI(idx);
            try
                save(roiFile,'roi_corners','-v7.3');
            catch
                warning('Could not save ROI file:');
            end

            % Return focus for more keys
            try
                uicontrol(focusCtl);
            catch
                figure(f);
            end

        %}
        case 'm'  % modify (start an editable rectangle initialized to this ROI)
            startEdit(idx);

        otherwise
            warndlg('Action must be "n" (remove) or "m" (modify).','ROI Action','modal');
    end
end

function startEdit(idx)
    % Start an editor rectangle at an existing ROI's position
    if idx < 1 || idx > size(roi_corners,1)
        return
    end

    % If an editor is already active, cancel it first
    if ~isempty(curRect) && isgraphics(curRect)
        delete(curRect);
        curRect = [];
    end

    % Seed the editor from current stored corners
    [x,y] = cornersToXY(roi_corners{idx,1});
    [pos, ~, ~] = cornersToRect([x,y]);

    % Enable minimal interactions for ROI editing
    try
        toggleTracking(f,ax,true);
    catch
    end

    curRect = drawrectangle(ax, ...
        'Position', pos, ...
        'Color', 'y', ...
        'LineWidth', 3, ...
        'FaceAlpha', 0.05);
    activeStarted = true;
    editTargetIdx = idx;         % <-- tell commit to REPLACE this index

    % Bring keyboard focus back to this figure
    try
        uicontrol(focusCtl);
    catch
        figure(f);
    end
    updateTitle();
end

function removeROI(idx)
    % Delete graphics for this ROI
    if idx <= numel(rectH) && isvalidHandle(rectH(idx))
        delete(rectH(idx));
        rectH(idx) = gobjects(1);
    end
    if idx <= numel(labelH) && isvalidHandle(labelH(idx))
        delete(labelH(idx));
        labelH(idx) = gobjects(1);
    end

    % Remove data
    if ~isempty(roi_corners)
        roi_corners(idx,:) = [];
    end

    % Shift graphics arrays down to stay aligned with roi_corners
    rectH(idx:end-1)  = rectH(idx+1:end);
    labelH(idx:end-1) = labelH(idx+1:end);
    if ~isempty(rectH), rectH(end) = gobjects(1); end
    if ~isempty(labelH), labelH(end) = gobjects(1); end

    % Renumber labels and reposition centers based on stored corners
    for i = 1:numel(rectH)
        if isvalidHandle(rectH(i))
            % Update label position based on current rect position
            pos = rectH(i).Position;
            [cx,cy] = rectCenter(pos);

            if isvalidHandle(labelH(i))
                set(labelH(i),'String',sprintf('%03d',i),'Position',[cx,cy,0]);
            else
                labelH(i,1) = text(ax, cx, cy, sprintf('%03d',i), 'Color','c','FontSize',14, ...
                                   'HorizontalAlignment','center','FontName','Trebuchet MS','FontWeight','bold');
            end
        end
    end

    % Update next index
    nextIndex = size(roi_corners,1) + 1;

    updateTitle();
end


function deleteAndRenumberCellFolders(idx, ROIs_dir)
    % Delete "Cell_idx", then shift Cell_(k) -> Cell_(k-1) for k > idx.
    % Robust to missing gaps; only acts on folders named Cell_\d\d\d.
    cellDirs = dir(fullfile(ROIs_dir,'Cell_*'));
    cellDirs = cellDirs([cellDirs.isdir]);

    % Parse numbers
    nums = [];
    for ii = 1:numel(cellDirs)
        t = regexp(cellDirs(ii).name, '^Cell_(\d{3})$', 'tokens', 'once');
        if ~isempty(t)
            nums(end+1) = str2double(t{1}); %#ok<AGROW>
        end
    end
    if isempty(nums)
        return
    end
    nums = sort(nums);  % ascending

    % 1) Delete the target folder if it exists
    tgtName = sprintf('Cell_%03d', idx);
    tgtPath = fullfile(ROIs_dir, tgtName);
    if isfolder(tgtPath)
        % remove contents recursively
        rmdir(tgtPath, 's');
    end

    % 2) Shift higher numbers down by 1: Cell_(k) -> Cell_(k-1) for k > idx
    %    Do in ascending order (Cell_(idx+1) first), which is safe because
    %    Cell_(k-1) was just freed by previous step.
    higher = nums(nums > idx);
    for k = higher
        oldName = sprintf('Cell_%03d', k);
        newName = sprintf('Cell_%03d', k-1);
        oldPath = fullfile(ROIs_dir, oldName);
        newPath = fullfile(ROIs_dir, newName);

        movefile(oldPath, newPath);

        %{
        % These checks cannot logically happen - KLS 20250927
        if isfolder(oldPath)
            % If destination unexpectedly exists, try to remove it (rare)
            if isfolder(newPath)
                rmdir(newPath,'s');
            end
            movefile(oldPath, newPath);
        end
        %}
    end
end

function finalizeAndExport()
    fig2 = figure('Color','k');
    ax2 = axes('Parent',fig2);
    axis(ax2,'image');
    ax2.Visible = 'off';
    imshow(imgMin, IRM_LUT, 'Parent',ax2, 'InitialMagnification','fit');
    hold(ax2,'on');
    for ii = 1:size(roi_corners,1)
        [x,y] = cornersToXY(roi_corners{ii,1});
        [x,y] = stayWithinImage(size(imgMin), x, y);
        [pos, cx, cy] = cornersToRect([x,y]);
        rectangle(ax2, 'Position',pos, 'EdgeColor','y', 'LineWidth',2);
        tcol = 'c';
        if ii == size(roi_corners,1)
            tcol = 'r';
        end
        text(ax2, cx, cy, sprintf('%03d',ii), 'Color',tcol,'FontSize',18, ...
             'HorizontalAlignment','center','FontName','Trebuchet MS');
    end
    hold(ax2,'off');
    
    try
        exportgraphics(fig2, fullfile(Save_individual_acq_dir,'Manually_Selected_ROIs.png'), 'Resolution',300);
    catch
        warning('Could not export ROI overview: %s');
    end
    
    if isgraphics(fig2)
        close(fig2);
    end
end
end

% ===== file-local helpers (no shared workspace) =====
function [xr,yr,xb,yb] = collectBoundaryPolylines(B, rows, cols)
    xr = [];
    yr = [];
    xb = [];
    yb = [];
    for k = 1:numel(B)
        boundary = B{k};
        touches = any(boundary(:,1)==1 | boundary(:,1)==rows | boundary(:,2)==1 | boundary(:,2)==cols);
        if touches
            xr = [xr; boundary(:,2); NaN]; %#ok<AGROW>
            yr = [yr; boundary(:,1); NaN];
        else
            xb = [xb; boundary(:,2); NaN];
            yb = [yb; boundary(:,1); NaN];
        end
    end
end

function [cx,cy] = rectCenter(pos)
    % pos = [x y w h]
    cx = pos(1) + pos(3)/2;
    cy = pos(2) + pos(4)/2;
end

function tf = isvalidHandle(h)
    tf = ~isempty(h) && isgraphics(h);
end

function [posClamped, xyCorners] = clampRectToImage(pos, sz)
    H = sz(1);
    W = sz(2);
    x1 = pos(1);
    y1 = pos(2);
    x2 = pos(1) + pos(3);
    y2 = pos(2) + pos(4);
    
    if x2 < x1
        t  = x1;
        x1 = x2;
        x2 = t;
    end
    if y2 < y1
        t  = y1;
        y1 = y2;
        y2 = t;
    end
    
    x1 = max(1, min(x1, W));
    x2 = max(1, min(x2, W));
    y1 = max(1, min(y1, H));
    y2 = max(1, min(y2, H));
    
    x1 = round(x1);
    x2 = round(x2);
    y1 = round(y1);
    y2 = round(y2);
    
    posClamped = [x1, y1, max(1,x2-x1), max(1,y2-y1)];
    xyCorners  = [x1,y1; x2,y2];
end

function [pos, cx, cy] = cornersToRect(xy)
    x1 = min(xy(:,1));
    y1 = min(xy(:,2));
    x2 = max(xy(:,1));
    y2 = max(xy(:,2));
    pos = [x1, y1, x2-x1, y2-y1];
    cx  = (x1 + x2)/2;
    cy  = (y1 + y2)/2;
end

function [x,y] = cornersToXY(xy)
    x = xy(:,1);
    y = xy(:,2);
end

function [x,y] = stayWithinImage(sz, x, y)
    H = sz(1);
    W = sz(2);
    x = max(1, min(round(x), W));
    y = max(1, min(round(y), H));
end

function toggleTracking(fig, ax, onoff)
    if onoff   % enable minimal interactions (for ROI drawing)
        enableDefaultInteractivity(ax);
        try, iptPointerManager(fig,'enable'); end
        fig.WindowButtonMotionFcn = [];  % let ROI manage motion
    else       % disable everything (idle waiting for keys)
        disableDefaultInteractivity(ax);
        datacursormode(fig,'off');
        fig.WindowButtonMotionFcn = [];
        fig.WindowScrollWheelFcn  = [];
        fig.WindowButtonDownFcn   = [];
        try, iptPointerManager(fig,'disable'); end
    end
end