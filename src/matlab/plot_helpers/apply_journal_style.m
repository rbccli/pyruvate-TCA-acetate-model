function apply_journal_style(ax, cb, varargin)
% APPLY_JOURNAL_STYLE  Apply a unified journal plotting style to axes/colorbar
%
% Usage:
%   apply_journal_style(ax, cb)
%   apply_journal_style(ax, [])
%   apply_journal_style(ax, cb, 'FontName', 'Helvetica', 'AxisFontSize', 10)
%
% Inputs:
%   ax   - axes handle
%   cb   - colorbar handle, or [] if no colorbar
%
% Optional name-value pairs:
%   'FontName'         - default 'Arial'
%   'AxisFontSize'     - default 9
%   'LabelFontSize'    - default 10
%   'CBFontSize'       - default 9
%   'CBLabelFontSize'  - default 10
%   'LineWidth'        - default 1
%   'TickDir'          - default 'out'
%   'Box'              - default 'off'
%   'TitleFontWeight'  - default 'normal'

    if nargin < 2
        cb = [];
    end

    p = inputParser;
    addParameter(p, 'FontName', 'Arial');
    addParameter(p, 'AxisFontSize', 9);
    addParameter(p, 'LabelFontSize', 10);
    addParameter(p, 'CBFontSize', 9);
    addParameter(p, 'CBLabelFontSize', 10);
    addParameter(p, 'LineWidth', 1);
    addParameter(p, 'TickDir', 'out');
    addParameter(p, 'Box', 'off');
    addParameter(p, 'TitleFontWeight', 'normal');
    parse(p, varargin{:});
    s = p.Results;

    % ---- axes style ----
    set(ax, ...
        'FontName', s.FontName, ...
        'FontSize', s.AxisFontSize, ...
        'LineWidth', s.LineWidth, ...
        'TickDir', s.TickDir, ...
        'Box', s.Box, ...
        'Layer', 'top');

    % ---- axis labels ----
    if isgraphics(ax.XLabel)
        ax.XLabel.FontName = s.FontName;
        ax.XLabel.FontSize = s.LabelFontSize;
    end

    if isgraphics(ax.YLabel)
        ax.YLabel.FontName = s.FontName;
        ax.YLabel.FontSize = s.LabelFontSize;
    end

    % ---- title ----
    if isgraphics(ax.Title)
        ax.Title.FontName = s.FontName;
        ax.Title.FontSize = s.LabelFontSize;
        ax.Title.FontWeight = s.TitleFontWeight;
    end

    % ---- colorbar style ----
    if ~isempty(cb) && isgraphics(cb)
        set(cb, ...
            'FontName', s.FontName, ...
            'FontSize', s.CBFontSize, ...
            'LineWidth', s.LineWidth, ...
            'TickDirection', s.TickDir);

        if isgraphics(cb.Label)
            cb.Label.FontName = s.FontName;
            cb.Label.FontSize = s.CBLabelFontSize;
        end
    end
end