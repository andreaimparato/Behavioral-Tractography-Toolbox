function [R2,P2]=MIXED_CORR_MATRIX(DATA,SUBJECT_ID) 
% This function computes the correlation matrix between columns of the input 
% data matrix 'DATA', while considering subject IDs provided in 'SUBJECT_ID'.
% 
% INPUTS:
%   - DATA: A matrix where each column represents a different variable or 
%           measurement, and each row represents an observation.
%   - SUBJECT_ID: A vector containing subject IDs corresponding to the rows
%                 of 'DATA'.
%
% OUTPUTS:
%   - R2: Correlation matrix between variables.
%   - P2: P-values corresponding to the correlations in the matrix.

clear R2
clear P2

% Loop through each column of the data matrix
for iter=1:size(DATA,2)
    iter
    % Loop through each column of the data matrix again
    for iter2=1:size(DATA,2)
        
        % Standardize the first measurement column
        MEAS=(DATA(:,iter));
        X=MEAS;
        if any(isnan(X(:)))
            xmu=nanmean(X);
            xsigma=nanstd(X);
            MEAS=(X-repmat(xmu,length(X),1))./repmat(xsigma,length(X),1);
        else
            [MEAS,xmu,xsigma]=zscore(X);
        end

        % Standardize the second measurement column
        MEAS2=(DATA(:,iter2));
        X=MEAS2;
        if any(isnan(X(:)))
            xmu=nanmean(X);
            xsigma=nanstd(X);
            MEAS2=(X-repmat(xmu,length(X),1))./repmat(xsigma,length(X),1);
        else
            [MEAS2,xmu,xsigma]=zscore(X);
        end
        
        % Retrieve subject IDs
        ID_F=SUBJECT_ID;

        % Create a dataset for linear mixed effects model
        ds=dataset(MEAS,MEAS2,ID_F);
        
        % Fit a linear mixed effects model
        % Commenting out this line as it's not being used currently
%         lme=fitlme(ds, 'MEAS~MEAS2 +(1|ID_F)+(MEAS2-1|ID_F)');
        
        % Fit a linear mixed effects model considering only fixed effects
        lme=fitlme(ds, 'MEAS~MEAS2 +(1|ID_F)');
        
        % Store the coefficient corresponding to the second variable as R2
        R2(iter,iter2)=double(lme.Coefficients(2,2));
        
        % Store the p-value corresponding to the coefficient as P2
        P2(iter,iter2)=double(lme.Coefficients(2,6));

    end
end
