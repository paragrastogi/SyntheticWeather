function [TSB,TSBidx] = BootstrapCustom(TS, B, nboot, SubPeriods)
% % Bootstrap function First check if the incoming time
% series can be split into B-sized blocks

if isrow(TS)
    TS = TS';
end

N = length(TS);

if mod(length(TS),B)>0
    % Number of blocks
    noB = floor(length(TS)/B)+1;
    % 'Missing' length that needs to be added to the
    % original series to make a complete end block
    AugLen = B - mod(length(TS),B);
    % Augment time series
    TSaug = [TS; TS(1:AugLen)];
    % Augmented index
    Idxaug = [(1:N)'; (1:AugLen)'];
else
    noB = (length(TS)/B);
    AugLen = 0;
    TSaug = TS;
    Idxaug = (1:N)';
end

fprintf(['Number of blocks for bootstrapping are ', ...
    '%d.\r\n The incoming time series was augmented', ...
	' with %d values\r\n'], noB, AugLen);

% Form the blocks from the augmented TS
Blocks = reshape(TSaug,  B, []);
% Store the original indices
OrigIdx = reshape(Idxaug,  B, []);

% Check if splitting the data into blocks leaves a remainder
% or not
if mod(size(Blocks,2),SubPeriods)>0
    % Size of the sub-block
    SubBloSi = floor(size(Blocks,2)/SubPeriods);
    % The left-over
    SubBloRem = mod(size(Blocks,2),SubPeriods);
else
    SubBloSi = (size(Blocks,2)/SubPeriods);
    SubBloRem = 0;
end

% Initialise the bootstrapped values and indices
TSBidx = NaN(size(TS,1),nboot); TSB = NaN(size(TS,1),nboot);

for m = 1:SubPeriods
    
    % Choose the current sub-period
    SubPerPicker = ((m-1)*SubBloSi)+1:( ...
		(m-1)*SubBloSi)+SubBloSi;
    CurrSubPer = Blocks(:, SubPerPicker);
    % These are its corresponding indices
    CurrIdx = OrigIdx(:, SubPerPicker);
    
    % Now bootstrap the blocks in this sub-period
    for k = 1:nboot
        
        % % Sample the data with replacement
        [Resampled, BSidx] = datasample(CurrSubPer, ...
            size(CurrSubPer,2), 2, 'Replace',true);
        % Save the Indices from the resample run as well.
        ResampledIdx = CurrIdx(:,BSidx);
        
        % Reshape the resampled blocks into a vector
        Resampled_Reshaped = reshape(Resampled,[],1);
        ResampledIdx_Reshaped = reshape(ResampledIdx,[],1);
        LenRes = length(Resampled_Reshaped);
        
        TSB(((m-1)*LenRes)+1:((m-1)*LenRes)+LenRes,k) = ...
            Resampled_Reshaped;
        TSBidx(((m-1)*LenRes)+1:((m-1)*LenRes)+LenRes,k) = ...
            ResampledIdx_Reshaped;
        
        if m == SubPeriods
            
            % Use the first hours of the series to fill in
            % the missing hours at the end of the time
            % series
            
            LenMiss = sum(isnan(TSB(:,k)));
            
            TSB(end-LenMiss+1:end,k) = TSB(1:LenMiss,k);
            TSBidx(end-LenMiss+1:end,k) = TSBidx(1:LenMiss,k);
            
        end
        
    end
    
end

TSB = TSB(1:length(TS),:);
TSBidx = TSBidx(1:length(TS),:);

end