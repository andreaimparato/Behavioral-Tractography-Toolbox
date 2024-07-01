function [R2, P2] = cross_corr(Symptoms, F, ParticipantID)
% Function to compute cross-correlation between symptoms for different participants
% Inputs:
%   Symptoms: Matrix of symptom data, where each column represents a symptom and each row represents a participant
%   F: Indices of the participants for which the correlation is to be computed
%   ParticipantID: Participant IDs
% Outputs:
%   R2: Matrix of correlation coefficients
%   P2: Matrix of p-values

for iter=1:size(Symptoms,2)
    iter; % Display the current iteration (not affecting the functionality)
    for iter2=1:size(Symptoms,2)
        
        % Standardize the first symptom data
        MEAS=(Symptoms(F,iter));
        X=MEAS;
        if any(isnan(X(:)))
            xmu=nanmean(X);
            xsigma=nanstd(X);
            MEAS=(X-repmat(xmu,length(X),1))./repmat(xsigma,length(X),1);
        else
            [MEAS,xmu,xsigma]=zscore(X);
        end

        % Standardize the second symptom data
        MEAS2=(Symptoms(F,iter2));
        X=MEAS2;
        if any(isnan(X(:)))
            xmu=nanmean(X);
            xsigma=nanstd(X);
            MEAS2=(X-repmat(xmu,length(X),1))./repmat(xsigma,length(X),1);
        else
            [MEAS2,xmu,xsigma]=zscore(X);
        end
        
        ID_F=ParticipantID(F,:);
        
        % Create a dataset with standardized symptom data and participant IDs
        ds=dataset(MEAS,MEAS2,ID_F);
        
        % Fit a linear mixed-effects model
        lme=fitlme(ds, 'MEAS~MEAS2 +(1|ID_F)');
        
        % Store the correlation coefficient and p-value
        R2(iter,iter2)=double(lme.Coefficients(2,2));
        P2(iter,iter2)=double(lme.Coefficients(2,6));
    end
end
