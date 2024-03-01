function GenerateTransportFiles

cd(getenv("froot_tools"));

% Define array of experiment IDs for which to generate transport files
runID = ["PTDC_003"];

% Specify whether the transport is calculated for the ice front or an
% ice-shelf draft contour line. The length of this array is the same as
% runID
location = ["icefront"]; % values are "icefront" or "draft"

% If location = draft, specify if the contour is moving with the evolving
% draft, or the location remains fixed at the original contour line
moving = 1; % can be 1 or 0

% Specify the depth above/below which properties should be calculated. If
% location = draft, this will also define the depth of the draft contour.
% The length of this array is the same as runID
depth = -400;


for ii=1:numel(runID)

    if location(ii) == "icefront"
        figurefiletag = location+"_below"+string(abs(depth(ii)))+"m";
    elseif location(ii) == "draft"
        if moving(ii)
            figurefiletag = "moving"+string(abs(depth(ii)))+"mdraft";
        else
            figurefiletag = "fixed"+string(abs(depth(ii)))+"mdraft";
        end
    else
        error("unknown location");
    end

     fprintf("Generating heatvolumetransport_%s_%s\n",figurefiletag,runID(ii));
     TSUVQMsection(depth(jj),figurefiletag(jj),runID(ii),0);
end
