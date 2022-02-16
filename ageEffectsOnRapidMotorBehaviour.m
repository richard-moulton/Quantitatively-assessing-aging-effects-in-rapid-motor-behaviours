%% Quantitatively assessing aging effects in rapid motor behaviours
%  Analyses data and generates figures for the paper "Quantitatively 
%  assessing aging effects in rapid motor behaviours: a cross-sectional
%  study" Moulton et al., 2022. Submitted to the Journal of
%  NeuroEngineering and Rehabilitation.
%  Copyright (C) 2022  Richard Hugh Moulton
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%  This script proceeds in five sections.
%   - First, global flags are set to control the script's output.
%   - Second, the data is loaded and variable names are initialized.
%   - Third, the Object-Hit data is analysed.
%   - Fourth, the Object-Hit-and-Avoid data is analysed.
%   - Finally, supporting functions appear at the end of the script.
%
%  Dependencies:
%   - savefigas.m (included in this directory)
%   - cohensD.m (included in this directory)
%   - The Statistics and Machine Learning Toolbox
%
%  14 February, 2022

%#ok<*ST2NM>

%% Global Flags
verbose = 1;            % Print to console throughout execution
graphFlag = 1;          % Produce figures throughout
fileSaveType = 'png';   % File format for all produced figures

%% Load Data and Initialize Variable Names
ohControl = load('Quantitatively-assessing-aging(OH).mat');
ohaControl = load('Quantitatively-assessing-aging(OHA).mat');
demographics = load('Quantitatively-assessing-aging_Participant-demographics.mat');

varNames = {'Steady-state rate','Targets hit','Mean hand speed left','Mean hand speed right','Mean hand speed bias',...
    'Movement area left','Movement area right','Movement area bias','Hand bias of hits','Hand selection overlap',...
    'Hand transition','Miss bias','Median error','Distractor proportion','Object processing rate','Task score'};

%% Object-Hit
% Build and export the necessary data
ohData = buildData(ohControl,demographics,varNames);
fileNameString = sprintf('ageSexEffects-OH-SSR.csv');
writematrix(ohData,fileNameString);

% Perform age-based demographic breakdown of participants
demographicAnalysis(ohData,'OH');

if graphFlag
    % See the range of targets hit achieved by participants
    f = rangeOfHitsHistogram(ohControl.hitsScores,'Object-Hit');
    fileNameString = sprintf('ageSexEffects-OH-RangeOfHits');
    savefigas(f,fileNameString,fileSaveType);
    
    % See how accurate participants are by horizontal bin during their early
    % and overwhelmed phases
    f = binAccuracyPlots(ohControl.earlyBinAccuracies,ohControl.overwhelmedBinAccuracies);
    fileNameString = sprintf('ageSexEffects-OH-PerBinAccuracies');
    savefigas(f,fileNameString,fileSaveType);
end

% Perform linear regressions between age and each variable
regressions = performAgeRegressions(ohData,varNames,'OH',fileSaveType);
fileNameString = sprintf('ageSexEffects-OH-AgeRegressions.mat');
save(fileNameString,'regressions');

% Plot the linear regressions between targets hit and steady-state rates
if graphFlag
    [f,mdl] = scatterWithRegression(ohControl.steadyStateRates,ohControl.hitsScores,[0 4],[100 300],0,1);
    fileNameString = sprintf('ageSexEffects-OH-SSRvHits');
    savefigas(f,fileNameString,fileSaveType);
    if verbose
        fprintf('OH Regression Model: Target Hits = %0.4f * Steady-state rate + %0.4f\n\n',...
            mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
    end
end

% Analyse whether any variable shows a significant difference between male
% and female participants
analysis = performSexAnalysis(ohData,varNames,'OH',verbose);
fileNameString = sprintf('ageSexEffects-OH-SexTtests.mat');
save(fileNameString,'analysis');

fprintf('===========================================================================================\n\n');

%% Object-Hit-and-Avoid
% Build and export the necessary data
ohaData = buildData(ohaControl,demographics,varNames);
fileNameString = sprintf('ageSexEffects-OHA-SSR.csv');
writematrix(ohaData,fileNameString);

% Perform age-based demographic breakdown of participants
demographicAnalysis(ohaData,'OHA');

if graphFlag
    % See the range of targets hit achieved by participants
    f = rangeOfHitsHistogram(ohaControl.hitsScores,'Object-Hit-and-Avoid');
    fileNameString = sprintf('ageSexEffects-OHA-RangeOfHits');
    savefigas(f,fileNameString,fileSaveType);
    
    % See how accurate participants are by horizontal bin during their early
    % and overwhelmed phases
    f = binAccuracyPlots(ohaControl.earlyBinAccuracies,ohaControl.overwhelmedBinAccuracies);
    fileNameString = sprintf('ageSexEffects-OHA-PerBinAccuracies');
    savefigas(f,fileNameString,fileSaveType);
    
    % Investigate statistics relating to the introduction of distractors
    [f,mdl] = scatterWithRegression(ohaControl.hitsScores/200,ohaControl.dropsScores/100,[0 1], [0 1], 3, 1);
    ci = coefCI(mdl);
    fprintf('For all subjects, the regression between hit proportion and drop proportion has a slope of %0.4f (%0.4f, %0.4f), intercept of %0.4f, and R2 of %0.4f.\n',...
        mdl.Coefficients.Estimate(2),ci(2,1),ci(2,2),mdl.Coefficients.Estimate(1),mdl.Rsquared.Adjusted);
    fileNameString = sprintf('ageSexEffects-OHA-HitsDrops(All)');
    savefigas(f,fileNameString,fileSaveType);
end

% Determine whether male and female participants respond to distractor
% objects differently
femaleIdx = find(ohaData(:,3)==0);
maleIdx = find(ohaData(:,3)==1);

[f,mdl] = scatterWithRegression(ohaControl.hitsScores(femaleIdx),ohaControl.dropsScores(femaleIdx),[0 200], [0 100], 2, 1);
ci = coefCI(mdl);
fprintf('For female subjects, the regression between hits and drops has a slope of %0.4f (%0.4f, %0.4f).\n',...
    mdl.Coefficients.Estimate(2),ci(2,1),ci(2,2));
fileNameString = sprintf('ageSexEffects-OHA-HitsDrops(F)');
savefigas(f,fileNameString,fileSaveType);

[f,mdl] = scatterWithRegression(ohaControl.hitsScores(maleIdx),ohaControl.dropsScores(maleIdx),[0 200], [0 100], 2, 1);
ci = coefCI(mdl);
fprintf('For male subjects, the regression between hits and drops has a slope of %0.4f (%0.4f, %0.4f).\n',...
    mdl.Coefficients.Estimate(2),ci(2,1),ci(2,2));
fileNameString = sprintf('ageSexEffects-OHA-HitsDrops(M)');
savefigas(f,fileNameString,fileSaveType);

% Perform linear regressions between age and each variable
regressions = performAgeRegressions(ohaData,varNames,'OHA',fileSaveType);
fileNameString = sprintf('ageSexEffects-OHA-AgeRegressions.mat');
save(fileNameString,'regressions');

% Plot the linear regressions between targets hit and steady-state rates
if graphFlag
    [f,mdl] = scatterWithRegression(ohaControl.steadyStateRates,ohaControl.hitsScores,[0 2],[50 200],0,1);
    fileNameString = sprintf('ageSexEffects-OHA-SSRvHits');
    savefigas(f,fileNameString,fileSaveType);
    fprintf('OHA Regression Model: Target Hits = %0.4f * Steady-state rate + %0.4f\n\n',...
        mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
end

analysis = performSexAnalysis(ohaData,varNames,'OHA',verbose);
fileNameString = sprintf('ageSexEffects-OHA-SexTtests.mat');
save(fileNameString,'analysis');

%% FUNCTIONS
% Build a single matrix containing all of the trials as rows and all of the
% variables of interest as columns.
function data = buildData(trialData,demographics,varNames)
numTrials = length(trialData.trialNames);
data = NaN(numTrials,length(varNames)+4);

for i=1:numTrials
    subjectID = trialData.subjectIDs(i);
    subject = demographics.subjects(subjectID);
    switch subject{1}
        case 'F'
            subjectSex = 0;
        case 'f'
            subjectSex = 0;
        case 'M'
            subjectSex = 1;
        case 'm'
            subjectSex = 1;
        otherwise
            subjectSex = 2;
    end
    % Build ID, Age, Sex, Handedness, SSR, Metadata, Hits
    subjectData = [subjectID getAge(trialData.trialNames(i),subject{2}) subjectSex str2num(subject{3}) trialData.steadyStateRates(i) trialData.hitsScores(i)];
    data(i,:) = [subjectData getMetadata(trialData.metaDataArray(i),varNames(3:end-1)) trialData.taskScores(i)];
end
end

% Perform linear regressions between age and various KINARM-measured
% variables. Save the resulting scatter-and-regression figures as well as
% the regression statistics.
function regressions = performAgeRegressions(data,varNames,taskName,fileSaveType)
regressions = NaN(length(varNames),13);
xLimits = [floor(min(data(:,2))/10)*10 ceil(max(data(:,2)/10))*10];
yLimits = [0 4; 100 300; 0 1; 0 1; -0.5 0.5; 0 0.25; 0 0.25; -0.5 0.5; -0.5 0.5; 0 0.4; -0.1 0.1; -0.25 0.25; 50 90; 0 35; 1 3.5; 0 5];
if strcmp(taskName,'OHA')
    yLimits(1,:) = [0 2];
    yLimits(2,:) = [0 200];
end

for i = 5:size(data,2)
    if ~isnan(data(1,i))
        varIdx = i-4;
        unityFlag = 0;
        ciFlag = 1;
        [f,mdl] = scatterWithRegression(data(:,2),data(:,i),xLimits,yLimits(varIdx,:),unityFlag,ciFlag);
        fileSaveName = sprintf('ageSexRegression-%s-Age-%s',taskName,varNames{varIdx});
        savefigas(f,fileSaveName,fileSaveType);
        
        ci = coefCI(mdl);
        
        regressions(varIdx,1) = mdl.Coefficients.Estimate(1);
        regressions(varIdx,2:3) = ci(1,:);
        regressions(varIdx,4) = mdl.Coefficients.pValue(1);
        regressions(varIdx,5) = mdl.Coefficients.tStat(1);
        regressions(varIdx,6) = mdl.Coefficients.Estimate(2);
        regressions(varIdx,7:8) = ci(2,:);
        regressions(varIdx,9) = mdl.Coefficients.pValue(2);
        regressions(varIdx,10) = mdl.Coefficients.tStat(2);
        regressions(varIdx,11) = mdl.RMSE;
        regressions(varIdx,12) = mdl.Rsquared.Adjusted;
        regressions(varIdx,13) = abs(10 * mdl.Coefficients.Estimate(2) / mean(data(:,i)));
    end
end

regressions = array2table(regressions,'RowNames',varNames,...
    'VariableNames',{'Intercept Estimate','Intercept 95% CI Low','Intercept 95% CI High','Intercept p-Value',...
    'Intercept t-Stat','Slope Estimate','Slope 95% CI Low','Slope 95% CI High','Slope p-Value','Slope t-Stat',...
    'Model RMSE','Model R^{2}','Decade Effect Size'});
end

% Perform t-tests to determine whether or not there are significant
% differences in parameters between male and female participants.
% Additionally, compute the Z statistic to determine whether the
% age-related changes in parameters are the same between male and female
% participants.
function analysis = performSexAnalysis(data,varNames,taskName,verbose)
maleData = data(data(:,3)==1,:);
femaleData = data(data(:,3)==0,:);
analysis = NaN(4,length(varNames));

for i=1:length(varNames)
    [~,p,ci] = ttest2(maleData(:,i+4),femaleData(:,i+4),'VarType','unequal');
    analysis(1,i) = p;
    analysis(2:3,i) = ci;
    analysis(4,i) = cohensD(maleData(:,i+4),femaleData(:,i+4));
end

analysis = array2table(analysis,'VariableNames',varNames,'RowNames',{'t-test p-value',...
    '95% CI Lower Bound','95% CI Upper Bound','Cohen''s d'});

% Check variables with significant age-related effects for a sex-based
% effect on these age-related effects.
sigVars = [5, 6, 17, 18, 19];
for i = 1:length(sigVars)
    idx = sigVars(i);
    if isnan(maleData(1,idx))
        continue
    end
    
    switch taskName
        case 'OH'
            switch idx
                case 5
                    yLimits = [0 4];
                case 6
                    yLimits = [100 300];
                case 17
                    yLimits = [50 100];
            end
        case 'OHA'
            switch idx
                case 5
                    yLimits = [0 2];
                case 6
                    yLimits = [50 200];
                case 17
                    yLimits = [50 90];
                case 18
                    yLimits = [0 35];
                case 19
                    yLimits = [1 3.5];
            end
    end
    
    [~,mdl] = scatterWithRegression(maleData(:,2),maleData(:,idx),[10 100],yLimits,0,1);
    bMaleEst = mdl.Coefficients.Estimate(2);
    bMaleSE = mdl.Coefficients.SE(2);
    if verbose
        fprintf('The %s v Age regression for male participants in %s has a slope of %0.4f with a standard error of %0.4f.\n',...
            varNames{idx-4},taskName,bMaleEst,bMaleSE);
    end
    
    [~,mdl] = scatterWithRegression(femaleData(:,2),femaleData(:,idx),[10 100],yLimits,0,1);
    bFemaleEst = mdl.Coefficients.Estimate(2);
    bFemaleSE = mdl.Coefficients.SE(2);
    if verbose
        fprintf('The %s v Age regression for female participants in %s has a slope of %0.4f with a standard error of %0.4f.\n',...
            varNames{idx-4},taskName,bFemaleEst,bFemaleSE);
    end
    
    Z = (bMaleEst - bFemaleEst)/sqrt(bMaleSE^2 + bFemaleSE^2);
    fprintf('For %s, the Z-statistic for %s regression slopes being equal by sex is %0.4f for a p-value of %0.4f.\n',...
        taskName,varNames{idx-4},Z,2*normcdf(-abs(Z)));
end
end

% Calculate a participant's age on the day of a trial by comparing their
% date of birth with the date of the trial.
function age = getAge(trialName,DOB)
trialDate = regexp(trialName,'(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_','names');
trialDate = datetime(str2num(trialDate.year),str2num(trialDate.month),str2num(trialDate.day));

dobDate = regexp(DOB,'(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})','names');
dobDate = datetime(str2num(dobDate.year),str2num(dobDate.month),str2num(dobDate.day));

age = split(between(dobDate,trialDate),{'years'});
end

% Make sure that meta data is extracted in the correct order, in case
% something happens with how metadata gets stored in the future.
function metaData = getMetadata(inputArray,varNames)
metaData = length(varNames);
for i=1:length(varNames)
    switch varNames{i}
        case 'Mean hand speed left'
            metaData(i) = inputArray.MEAN_HAND_SPEED_LEFT;
        case 'Mean hand speed right'
            metaData(i) = inputArray.MEAN_HAND_SPEED_RIGHT;
        case 'Mean hand speed bias'
            metaData(i) = inputArray.MEAN_HAND_SPEED_BIAS;
        case 'Movement area left'
            metaData(i) = inputArray.MOVEMENT_AREA_LEFT;
        case 'Movement area right'
            metaData(i) = inputArray.MOVEMENT_AREA_RIGHT;
        case 'Movement area bias'
            metaData(i) = inputArray.MOVEMENT_AREA_BIAS;
        case 'Hand bias of hits'
            metaData(i) = inputArray.HAND_BIAS_OF_HITS;
        case 'Hand selection overlap'
            metaData(i) = inputArray.HAND_SELECTION_OVERLAP;
        case 'Hand transition'
            metaData(i) = inputArray.HAND_TRANSITION;
        case 'Miss bias'
            metaData(i) = inputArray.MISS_BIAS;
        case 'Distractor proportion'
            if isfield(inputArray,'DISTRACTOR_PROPORTION')
                metaData(i) = inputArray.DISTRACTOR_PROPORTION;
            else
                metaData(i) = NaN;
            end
        case 'Median error'
            metaData(i) = inputArray.MEDIAN_ERROR;
        case 'Object processing rate'
            if isfield(inputArray,'OBJECT_PROCESSING_RATE')
                metaData(i) = inputArray.OBJECT_PROCESSING_RATE;
            else
                metaData(i) = NaN;
            end
        otherwise
            fprintf('Unrecognized Variable Name! Got %s and it doesn''t match any of the variables names we know.\n',varNames{i});
            keyboard
    end
end
end

% Produce a scatter plot between two variables overlaid with a linear
% regression between these same variables. Include the unity line if
% unityFlag is set.
function [f,mdl] = scatterWithRegression(x,y,xLimits,yLimits,unityFlag,ciFlag)
f = figure;
hold on;

% Scatter
scatter(x,y,150,'black','filled');
pause(1);

% Unity Line
switch unityFlag
    case 1
        minCoord = min([xLimits(1)  yLimits(1)]);
        maxCoord = max([xLimits(2) yLimits(2)]);
        plot([minCoord maxCoord],[minCoord maxCoord],'k','LineWidth',6)
    case 2
        plot([0 200],[0 100],'k','LineWidth',6);
    case 3
        plot([0 1],[0 1],'k','LineWidth',6);
        plot([1 0],[0 1],'--k','LineWidth',6);
end

% Regression
mdl = fitlm(x,y,'linear');
xVals = linspace(xLimits(1),xLimits(2));
[mdlVals,mdlCIs] = predict(mdl,xVals','Alpha',0.05);

plot(xVals,mdlVals,'Color',[0.8500 0.3250 0.0980],'LineWidth',5);
if ciFlag
    plot(xVals,mdlCIs(:,1),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    plot(xVals,mdlCIs(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2);
end

xlim(xLimits);
ylim(yLimits);

ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 15;
f.WindowState = 'maximized';
axis square
hold off;
end

% Determine the number of participants and their breakdown by age, sex, and
% handedness for a given task.
function demographicAnalysis(data,taskName)
fprintf('Conducting demographic analysis for %s.\n',taskName);

uniqueParticipants = unique(data(:,1));
totalParticipants = length(uniqueParticipants);
participantDemos = NaN(totalParticipants,3);

for i=1:totalParticipants
    participantIdx = find(data(:,1)==uniqueParticipants(i),1);
    % Age, Sex, Handedness
    participantDemos(i,:) = data(participantIdx,2:4);
end

ageRanges = [0 100;18 29;30 39;40 49;50 59;60 69;70 79;80 100];

for i=1:size(ageRanges,1)
    ageParticipants = participantDemos((participantDemos(:,1) >= ageRanges(i,1)) & (participantDemos(:,1) <= ageRanges(i,2)),:);
    femaleParticipants = ageParticipants(ageParticipants(:,2)==0,:);
    maleParticipants = ageParticipants(ageParticipants(:,2)==1,:);
    
    fprintf('For the age range %i-%i there are %i participants with a median age of %i.\n',...
        ageRanges(i,:),size(ageParticipants,1),median(ageParticipants(:,1)));
    fprintf('There are %i female participants and %i male participants.\n',size(femaleParticipants,1),size(maleParticipants,1));
    fprintf('There are %i/%i left-handed, %i/%i mixed-handed, and %i/%i right-handed participants (all F/M).\n\n',...
        length(find(femaleParticipants(:,3)==-10)),length(find(maleParticipants(:,3)==-10)),...
        length(find(femaleParticipants(:,3)==0)),length(find(maleParticipants(:,3)==0)),...
        length(find(femaleParticipants(:,3)==10)),length(find(maleParticipants(:,3)==10)));
end
end

% Produce a histogram of the first argument. Exclusively used for the
% number of targets hit in this script.
function f = rangeOfHitsHistogram(hitsScores,taskName)
% The type of task determines how many targets and participants are
% involved. The histogram's axis limits should be adjusted accordingly.
switch taskName
    case 'Object-Hit'
        maxHits = 300;
        maxSubj = 150;
    case 'Object-Hit-and-Avoid'
        maxHits = 200;
        maxSubj = 150;
    case 'Turbo Object-Hit'
        maxHits = 300;
        maxSubj = 700;
    case 'Turbo Object-Hit-and-Avoid'
        maxHits = 200;
        maxSubj = 700;
    otherwise
        fprintf('Unexpected task name! Got ''%s''.\n',taskName);
        keyboard
end

hitsRange = 0:(maxHits/20):maxHits;

f = figure;
hold on;
histogram(hitsScores,hitsRange,'FaceColor',[0 0 0],'FaceAlpha',0.6,'EdgeAlpha',0.6)
axis([0 maxHits 0 maxSubj]);
ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 10;
axis square
f.WindowState = 'maximized';
hold off;

% Normalize the hits scores and then perform the Kolmogorov-Smirnov test
% for normality
normHitsScores = (hitsScores-mean(hitsScores))/std(hitsScores);
[h,p] = kstest(normHitsScores);

if h
    fprintf('\nThe Kolmogorov-Smirnov test rejects the null hypothesis that the number of targets hit comes from the standard normal distribution (p = %0.4f).\n',p);
else
    fprintf('\nThe Kolmogorov-Smirnov test fails to reject the null hypothesis that the number of targets hit comes from the standard normal distribution (p = %0.4f).\n',p);
end

% Plot the empirical CDF and the standard normal CDF for comparison. 
cdfplot(normHitsScores)
hold on;
x_values = linspace(min(normHitsScores),max(normHitsScores));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')
hold off;

fprintf('\nThe number of targets hit has a range of %i-%i, a median value of %i, and an interquartile range of %0.2f.\n',...
    min(hitsScores),max(hitsScores),median(hitsScores),iqr(hitsScores));
end

% Produce scatter plots showing participants' average hitting accuracies by
% horizontal bin during the OH and OHA tasks. Can also be used to show
% participants' average distractor avoiding accuracies by horizontal bin.
function f = binAccuracyPlots(earlyBinAccuracies,overwhelmedBinAccuracies)

binAccuracies = NaN(10,2);
binSEM = NaN(10,2);
for i=1:size(binAccuracies,1)
    binAccuracies(i,1) = mean(earlyBinAccuracies(:,i),'omitNan');
    binSEM(i,1) = std(earlyBinAccuracies(:,i),'omitnan')/sqrt(sum(~isnan(earlyBinAccuracies(:,i))));
    binAccuracies(i,2) = mean(overwhelmedBinAccuracies(:,i),'omitNan');
    binSEM(i,2) = std(overwhelmedBinAccuracies(:,i),'omitnan')/sqrt(sum(~isnan(overwhelmedBinAccuracies(:,i))));
end

f = figure;
hold on;
colour1 = [0.0742 0.3867 0.6172];
colour2 = [0.6172 0.1445 0.0742];
errorbar(1:10,binAccuracies(:,1),binSEM(:,1),...
    'LineStyle','none','LineWidth',20,'CapSize',30,'Marker','o',...
    'MarkerSize',20,'MarkerFaceColor','auto','Color',[colour1 0.6]);
errorbar(1:10,binAccuracies(:,2),binSEM(:,2),...
    'LineStyle','none','LineWidth',20,'CapSize',30,'Marker','o',...
    'MarkerSize',20,'MarkerFaceColor','auto','Color',[colour2 0.6]);
xlim([0.5 10.5]);
ylim([0 1]);
xticks([1 4 7 10]);
yticks([0 1]);
ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 10;
axis square
f.WindowState = 'maximized';
hold off;
end
