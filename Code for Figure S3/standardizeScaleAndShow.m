function standardizeScaleAndShow(figs, saveFigs, figPath)
% this is a utility function to standardize the scale of groups of figures
% input is a cell array with each cell containing a group of at least two
% figures whose scales will be standardized relative to each other. this
% function also saves the figures to a specified path: figPath.
% By: Abdo Sharaf - abdo.sharaf@yale.edu. 
% Xinyuan Zheng - xinyuan.zheng@yale.edu. 
% Last Modified: July 23 2022
% *************************************************************************


% figs = figGroups
% figPath = svpath_ep
for grp = 1:length(figs) % grp = 3
    
    figGrp = figs{grp};
    
    if isempty(figGrp)
        continue;
    end
    
    if length(figGrp) < 2
        if length(findall(figGrp, 'Type', 'Axes')) < 2
            % set it visible 
            set(figGrp,'visible','on');
            
            if saveFigs
                % get its children axes
                fig_axs = findall(figGrp,'Type', 'Axes');
                % get and prep the save name for the figure
                savenm = get(fig_axs(1), 'Title');
                savenm = regexprep(savenm.String, ' ', '_');
                
                % save it
                
                % set size to full screen
%                 set(figGrp, 'Position', get(0, 'Screensize'));
                set(figGrp, 'Position', [100,100,700,500]);
                saveas(figGrp, fullfile(figPath, [savenm '.png']), 'png');
                saveas(figGrp, fullfile(figPath, [savenm '.eps']), 'eps');
                saveas(figGrp, fullfile(figPath, [savenm '.fig']), 'fig');
            end
                
            continue
        end
    end

    miny1 = 0;
    maxy2 = 0;
    allylims = []; 
    % loop once to find the standard y-axis limits
    for fig = 1:length(figGrp)
        fig_axs = findall(figGrp(fig),'Type', 'Axes');
        for ax = 1:length(fig_axs)
            ylims = ylim(fig_axs(ax));
            allylims = [allylims; ylims]; 
            if ylims(1) < miny1
                miny1 = ylims(1);
            end
            
            if ylims(2) > maxy2
                maxy2 = ylims(2);
            end
        end
    end
    
    allylims = sort(allylims);
    
    % now loop again to standardize
    for fig = 1:length(figGrp)
        fig_axs = findall(figGrp(fig),'Type', 'Axes');
        
        if isempty(fig_axs)
            continue;
        end
        
        % get and prep the save name for the figure 
        savenm = get(fig_axs(1), 'Title'); 
        savenm = regexprep(savenm.String, ' ', '_'); 
        
        % edit
        for ax = 1:length(fig_axs)
            ylims = ylim(fig_axs(ax));
            low = miny1; 
            high = maxy2; 
%             idx = size(allylims, 1); 
%             while high > ylims(2) + 500
%                 idx = idx - 1; 
%                 if idx <= 0
%                     break; 
%                 end
%                 high = allylims(idx, 2); 
%             end 
%             while low < ylims(1) - 500
%                 idx = idx - 1;
%                 if idx <= 0
%                     break;
%                 end
%                 low = allylims(idx, 2); 
%             end 
            
            try
                set(fig_axs(ax),'YLim',[low, high]);
            catch 
                warning(['A figure doesn''t have proper y-axis limits.' ...
                    ' Setting limits arbitrarily to 0 and 1']); 
                low = 0; high = 1; 
                set(fig_axs(ax), 'YLim', [low, high]); 
            end
        end
        
        if saveFigs
            % set size to full screen
            set(figGrp(fig), 'Position', [100,100,700,500]);
            saveas(figGrp(fig), fullfile(figPath, [savenm '.png']), 'png');
            saveas(figGrp(fig), fullfile(figPath, [savenm '.eps']), 'eps');
            saveas(figGrp(fig), fullfile(figPath, [savenm '.fig']), 'fig');
        end
        
        % show
        set(figGrp(fig),'visible','on');
    end
end
end