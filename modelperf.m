function [ stats, roc ] = modelperf( pred, target, T )
%MODELPERF Calculates many performance metrics used to determine the
% efficacy of a given set of predictions on a set of dichotomous or
% continuous targets.
%
%   [ STATS ] = modelperf(PRED,TARGET) calculates the following
%   measures and stores them in a field in the structure STATS:
%
%       AUROC   - Area under the receiver operator (characteristic) curve
%       HL      - Hosmer-Lemeshow Statistic
%       HLD     - Hosmer-Lemeshow Statistic normalized to variable range
%       acc     - Mean accuracy of the prediction using a threshold of 0.5
%       R       - Shapiro's R, the geometric mean of positive outcomes
%       B       - Brier's score, the mean square error
%       Badj    - Adjusted Brier's score, percent error of MSE to null MSE
%       SMR     - Standardized mortality ratio. Mean observed outcomes
%                   divided by mean predicted outcomes.
%       coxA    - Alpha coefficient in cox calibration.
%       coxB    - Beta coefficient in cox calibration.
%                   Cox calibration involves a linear regression of the
%                   predictions onto the targets. coxA is the constant
%                   offset, while coxB is the slope.
%       TP      - True positives.
%       FP      - False positives.
%       TN      - True negatives.
%       FN      - False negatives.
%       sens    - Sensitivity. Calculated as: TP/(TP+FN)
%       spec    - Specificity. Calculated as: TN/(TN+FP)
%       ppv     - Positive predictive value. Calculated as: TP/(TP+FP)
%       npv     - Negative predictive value. Calculated as: TN/(TN+FN)
%       ll      - Log likelihood
%       ber     - Balanced error rate
%
%   [ STATS, ROC ] = modelperf(PRED,TARGET) additionally calculates the
%       value pairs which form the ROC curve.
%           roc.x - Stores 1-specificity, plotted on the x-axis
%           roc.y - Stores Sensitivity, plotted on the y-axis
%           roc.t - Stores the threshold used for this value
%
%   [ STATS, ROC ] = modelperf(PRED,TARGET,T) uses the input threshold T in
%   the calculation of operating point statistics and measures derived from
%   operating point statistics (i.e. anything dependent on true positives,
%   true negatives, false positives, and false negatives).
%
%	Copyright 2015 Alistair Johnson

%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson
%	Contact: alistairewj@gmail.com

stats=[]; roc=[];

if nargin==0
    %TODO: Display this properly
    stats=idealStruct(0);
    fn = fieldnames(stats);
    fprintf('BINARY METRICS OUTPUT:\n');
    fprintf('Metric\tIdeal value\n');
    for f=1:numel(fn)
        fprintf('%s\t\t%g\n',fn{f},stats.(fn{f}));
    end
    stats=idealStruct(1);
    fn = fieldnames(stats);
    fprintf('CONTINUOUS METRICS OUTPUT:\n');
    fprintf('Metric\tIdeal value\n');
    for f=1:numel(fn)
        fprintf('%s\t\t%g\n',fn{f},stats.(fn{f}));
    end
    return;
end

if nargin<3
    T = 0.5;
end

%=== Add in some error checks on the input predictions
if min(size(pred))~=1 || min(size(target))~=1 || numel(pred)~=numel(target)
    error('modelperf:BadInput','This function expects vector inputs of the same size.');
end

% ensures the inputs are column vectors
pred = pred(:);
target = target(:);

if max(target>1) || max(target<0)
    CONTINUOUS_FLAG = 1;
else
    CONTINUOUS_FLAG = 0;
    if max(pred>1) || min(pred<0)
        warning('modelperf:IllConditioned',...
            ['\nThis function will return potentially erroneous values \n' ...
            'when predictions are not in the range [0,1].\n']);
    end
end

if isempty(pred)
    stats=idealStruct(CONTINUOUS_FLAG);
    fn = fieldnames(stats);
    for f=1:numel(fn)
        stats.(fn{f}) = [];
    end
    return;
elseif sum(isnan(pred))==length(pred)
    % Pred is all NaNs
    stats = idealStruct(CONTINUOUS_FLAG);
    fn = fieldnames(stats);
    for f=1:numel(fn)
        stats.(fn{f}) = [];
    end
    return;
end

%=========================================================%
%=== STATISTICS ON BOTH CONTINUOUS AND BINARY OUTCOMES ===%
%=========================================================%

% Shapiro's R (geometric mean of positive targetcomes)
% transform into naural logarithm, take arithmetic mean, transform back
stats.R=exp(sum(log(pred(target==1)))/sum(target));

% Brier's Score, B (mean square error)
stats.B=mean((target-pred).^2);
% Adjusted Brier's Score, also known as Efron's pseudo R squared or
% sum-of-squares R squared
stats.Badj = 1 - ( mean((target-pred).^2)  ) / ( mean( (target-mean(target)).^2 ) );
stats.r2efron = stats.Badj;

pred_binary = pred>=T;
stats.T = T;


if CONTINUOUS_FLAG == 1
    stats.r2 = corr(pred,target)^2;    
else
    
    % Accuracy
    % stats.acc = sum( (pred_binary >= T) == target ) / numel(target);
    stats.acc=1-(sum(abs(round(target-pred_binary)),1)/length(target));
    
    % hosmer-lemeshow
    stats.HL=lemeshow(pred,target);
    stats.HLD=lemeshowNormalized(pred,target);
    
    % Standardised Mortality Ratio (used in epidemiology)
    stats.SMR=sum(target)/sum(pred);
    
    idx1 = target==1;
    idx0 = target==0;
    stats.TP=sum(pred_binary(idx1));
    stats.FN=sum(~pred_binary(idx1));
    stats.FP=sum(pred_binary(idx0));
    stats.TN=sum(~pred_binary);
    
    stats.sens=stats.TP/(stats.TP+stats.FN);
    stats.spec=stats.TN/(stats.TN+stats.FP);
    stats.ppv=stats.TP/(stats.TP+stats.FP);
    stats.npv=stats.TN/(stats.TN+stats.FN);
    stats.ber=(stats.sens + stats.spec)/2;
    stats.f1 = 2*stats.sens*stats.ppv / (stats.sens+stats.ppv); % f1 score
    
    %=== Odds ratio!
    stats.OR = stats.TP*stats.TN / (stats.FN*stats.FP);
    
    stats.PositiveLR=(1-stats.sens)/stats.spec;
    stats.NegativeLR=(stats.sens)/(1-stats.spec);
    
    % negative log likelihood (NLL)
    stats.nll = sum(-target.*log(pred+eps)-(1-target).*log(1-pred+eps));
    
    % 1 - (NLL normalised to the null model) (i.e. McFadden's pseudo R2 or entropy R2)
    stats.r2entropy = 1 - ( sum(-target.*log(pred+eps)-(1-target).*log(1-pred+eps)) ) / sum(-target*log(mean(target)+eps)-(1-target).*log(1-mean(target)+eps));
    
    % NLL normalised to the number of data points
    stats.logloss = sum(-target.*log(pred+eps)-(1-target).*log(1-pred+eps)) / numel(target);
    
    %=== Calculate SeSp
    [th,idxSort] = sort(pred,1,'ascend');
    tar = target(idxSort);
    TP = flipud(tar);
    FP = cumsum(1-TP);
    FP = flipud(FP);
    TP = cumsum(TP);
    TP = flipud(TP);
    FN = cumsum(tar)-tar;
    TN = numel(tar) - TP - FP - FN;
    
    %=== Sensitivity (true positive rate)
    roc.y = TP ./ (TP + FN);
    %=== 1-Specificity (false positive rate)
    roc.x = 1- (TN ./ (TN + FP));
    
    idxSE = find(roc.y<0.5,1);
    idxSP = find(roc.x<0.5,1);
    stats.SeSp = 1-(roc.y(idxSP) - roc.x(idxSE) + 1)/2;
    
    % Minimum of the sensitivity and the PPV, defined as "s1"
    PPV = TP ./ (TP + FP);
    stats.s1 = min(stats.sens, stats.ppv);
    
    % Also calculate the best possible threshold for the "s1" measure
    [stats.s1best, stats.Tbest] =  max( min([PPV,roc.y],[],2), [], 1 );
    stats.Tbest = th(stats.Tbest);
    
    %============================%
    %=== Concordance measures ===%
    %============================%
    %=== Arrange predictions
    [predSorted,idx] = sort(pred,1,'ascend');
    targetSorted=target(idx);
    N = size(predSorted,1);
    
    %=== Find location of negative targets
    negative = targetSorted==0;
    
    %=== Count the number of negative targets below each element
    negativeCS = cumsum(negative,1);
    
    %=== Get number of negative elements below each positive outcome
    % Each element of pos represents the number of concordant pairs for that
    % positive outcome - this is used for the AUROC
    concordant = sum(negativeCS(~negative));
    
    %=== count number who are negative
    count = sum(negative,1);
    stats.AUROC = concordant ./ (count .* (N-count));
    
    % [ out.AUROCEXACT ] = wilcoxonEXACT( pred, target );
    
    % Gamma
    %=== Get number of positive elements above each positive outcome
    discordant = cumsum(targetSorted==1,1);
    discordant = sum(discordant(~negative));
    stats.gamma = (concordant - discordant) ./ (concordant+discordant);
    stats.gammaz = stats.gamma*sqrt((concordant+discordant)/(2*N*(1-stats.gamma^2))); % z-test, this stat is approx. normal
    
    % Find ties
    ties_x = diff(targetSorted)~=0 & diff(predSorted)==0;
    ties_y = diff(targetSorted)==0 & diff(predSorted)~=0;
    stats.dyx = (concordant - discordant) ./ (0.5*N*(N-1)-sum(ties_y));
    stats.dxy = (concordant - discordant) ./ (0.5*N*(N-1)-sum(ties_x));
    % % cox linear regression testing
    %     b=glmfit(pred,target,'binomial','link','identity');
    %     out.coxA=b(1); out.coxB=b(2);
end
end

function [ HL, H ] = lemeshow( pred, target, D )

if nargin<3
    D=10; % Number of deciles to use.
    % Note that the HL test is only defined for D=10 (or D=8 if applying to a training set)
end
if min(size(pred))~=1 || min(size(target))~=1
    error('Function only accepts vector inputs.');
end

% ensure they are column vectors
target = target(:); pred = pred(:);
N = numel(target);

[pred,idxSort] = sort(pred,1,'ascend');
target = target(idxSort);

% create an index of the outcomes into each of the D deciles
idxDecile = ceil( (1:N)/N*D )';

H = zeros(D,2);
H(:,1) = accumarray(idxDecile,target);
H(:,2) = accumarray(idxDecile,pred);

% calculate expected probability in each decile
n = hist(idxDecile,1:D);
p = H(:,2)./n;
HL = (H(:,1)-H(:,2)).^2 ./ (p.*n.*(1-p)+eps); % eps prevents division by 0
HL = sum(HL);


end


function [ HL ] = lemeshowNormalized( pred, target, D )

if nargin<3
    D=10;
end

% call lemeshow subfunction
[ HL, H ] = lemeshow( pred, target, D );

% normalise by the range of the predictions
if numel(H)>1
    HL = HL ./ (H(end)-H(1));
end

end

function [defStruct] = defaultStruct()
defStruct = struct('AUROC',[],...
    'HL',[],...
    'acc',[],...
    'SeSp',[],...
    'll',[],...
    'logloss',[],...
    'ber',[],...
    'f1',[],...
    'normll',[],...
    'T',[],...
    'R',[],...
    'B',[],...
    'SMR',[],...
    'coxA',[],...
    'coxB',[],...
    'TP',[],...
    'FP',[],...
    'TN',[],...
    'FN',[],...
    'sens',[],...
    'spec',[],...
    'ppv',[],...
    'npv',[],...
    'OR',[],...
    'PositiveLR',[],...
    'NegativeLR',[]);
end
function [defStruct] = idealStruct(CONTINUOUS_FLAG)
if CONTINUOUS_FLAG==1
    % display continuous struct
    defStruct = struct();
else % binary struct
    defStruct = struct('AUROC',1,...
        'HL',0,...
        'acc',1,...
        'SeSp',0,...
        'll',0,...
        'logloss',0,...
        'ber',0,...
        'f1',1,...
        'normll',0,...
        'T',[],...
        'R',0,...
        'B',0,...
        'SMR',1,...
        'coxA',0,...
        'coxB',1,...
        'TP',[],...
        'FP',[],...
        'TN',[],...
        'FN',[],...
        'sens',1,...
        'spec',1,...
        'ppv',1,...
        'npv',1,...
        'OR',[],...
        'PositiveLR',[],...
        'NegativeLR',[]);
end
end


function [ W ] = wilcoxonEXACT( pred, target )
%WILCOXON Calculates the Wilcoxon statistic, equivalent to the area under
%   the receiver operator characteristic (AUROC)
%
%   W = wilcoxon(pred, target) Calculates the Wilcoxon statistic W given a
%   vector of predictions ranging from 0-1 and their associates targets,
%   which are either 0 or 1
%
%	Inputs:
%       target - the target variable, 0 or 1, predicted by pred
%       pred   - the prediction made (the probability that target=1).
%   pred,target should be column vectors.
%
%	Outputs:
%		W - Probability( PRED|target=1 > PRED|target=0 )
%               Calculation: sum(sum(PRED|target=1 > PRED|target=0))
%               Equivalent to the AUROC.
%
%
%   Example 1: Classify ~90% correctly.
%       pred=rand(100,1);
%       target=round(pred);
%       target(1:10:end)=1-target(1:10:end);
%       W = wilcoxonEXACT(pred,target)

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2013-09-20 09:24:36 +0100 (Fri, 20 Sep 2013) $
%   $Revision: 850 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:11:42
%   Contact:alistairewj@gmail.com

% do not use nan targets in calculations
idxNotNan = isfinite(target);
idx=target==0;
negative=pred(idx & idxNotNan);
N1=length(negative);
positive=pred(~idx & idxNotNan);
N2=length(positive);

W=0; W2=0;
% compare 0s to 1s
for n=1:N2
    idx=gt(positive(n),negative);
    W=W+sum(idx);
    W2=W2+sum(eq(positive(n),negative)); % Double computation time, but more accurate.
end
W=W+W2*0.5;
W=W/(N1*N2);
end