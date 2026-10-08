function [fig, ax] = make_journal_figure(width_cm, height_cm)
% MAKE_JOURNAL_FIGURE  Create a white-background journal-style figure
%
% Usage:
%   [fig, ax] = make_journal_figure()
%   [fig, ax] = make_journal_figure(8.5, 7)
%
% Inputs:
%   width_cm   - figure width in cm (default 8.5)
%   height_cm  - figure height in cm (default 7)
%
% Outputs:
%   fig        - figure handle
%   ax         - axes handle

    if nargin < 1 || isempty(width_cm)
        width_cm = 8.5;
    end
    if nargin < 2 || isempty(height_cm)
        height_cm = 7;
    end

    fig = figure( ...
        'Color', 'w', ...
        'Units', 'centimeters', ...
        'Position', [5 5 width_cm height_cm]);

    ax = axes('Parent', fig);
end