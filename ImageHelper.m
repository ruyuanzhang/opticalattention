%{
    A class encapusalation of several methds useful for image alignment.
    All useful functions are static and expect trial structures / images as
    needed.

    Here, alignment refers to the temporal domain. This class is not
    concerned with motion correction or image registration.

    The goal is to create a 4D set of images X by Y by Time by Trial. 
%}


classdef ImageHelper
   
    properties(Constant)
        BAD_IMAGE_INT5 = [0 1 0 0 1]';
    end
    methods(Static) 
        %% Align All Physiologic Data
        % The new hotness - Align pupil, eye, and image data.
        function [E1,E2,P1,P2,I,ep_outtime,im_outtime] = resamplePhysiologicData(T,eyeImageTimeOffset,eyeSamplesPerImage)
            % UNFINISHED
            % Pupil and eye data
            [E1,E2] = ImageHelper.convertEyeData(T);
            P1 = double(T.pupilData(1:2:end));
            P2 = double(T.pupilData(2:2:end));
            
            % We assume that the first pupil/eye sample occurs at
            % TrialStart and the last occurs at TrialEnd. Looking at a
            % representative data file, it appears the average latency from
            % TrialStart to EyeData is 5 ms.
            e_intime = linspace(double(T.trialStartTime),double(T.trialEndTime),length(E1));
            p_intime = linspace(double(T.trialStartTime),double(T.trialEndTime),length(P1));
            
            pre = double(T.(time_zero_field)):-(100000/240):double(T.trialStartTime);
            post = double(T.(time_zero_field)):(100000/240):double(T.trialEndTime);
            
            ep_outtime = [fliplr(pre) post(2:end)];
            
            E1 = ImageHelper.interpolateAnything(E1,e_intime,ep_outtime);
            E2 = ImageHelper.interpolateAnything(E2,e_intime,ep_outtime);
            P1 = ImageHelper.interpolateAnything(P1,p_intime,ep_outtime);
            P2 = ImageHelper.interpolateAnything(P2,p_intime,ep_outtime);
            
            pre = double(T.(time_zero_field)):-(100000/5):double(T.trialStartTime);
            post = double(T.(time_zero_field)):(100000/5):(double(T.imageTime(end)));
            
        end
        
        function out = interpolateAnything(in, intime, outtime)
            % Wrapper for interpolations. Can be used to interplolate 1D
            % eye or pupil data, or imageData if properly shaped.
            interpdim = ndims(in);
            if(length(intime) ~= size(in,interpdim))
                error('Size mismatch between in and intime. Likely cause is a column instead of a row vector was provided.');
            end
            
            perm = 1:ndims(in);
            perm(1) = ndims(in);
            perm(end) = 1;
            
            in = permute(in,perm);
            out = interp1(intime,in,outtime);
            out = permute(out,perm);
        end
        
        %% Align images of a trial
        function [im,ix0] = parseImageStream(trials,images,timefieldname,contingency,parseLength)
            % Should behave almost identically to alignImages, except is
            % not agnostic as to what image belongs to what trial. Images
            % now only belong to the nearest instance of timefieldname
            % within +/- parseLength images
            if(nargin < 4)
                contingency = [];
            end
            if(nargin < 5)
                parseLength = 2.5 * 1000 * 100; % 5 seconds
            end
            
            FRAMETIME = 20000;
            numAlignmentPoints = length([trials.(timefieldname)]);
            numImages = 2*ceil(parseLength/FRAMETIME);
            ix0 = (numImages/2)+1;

            im = nan(size(images,1), size(images,2), numImages, numAlignmentPoints,'single');
            imdelete = false(1,1,1,size(im,4));
            
            eventTimes = ImageHelper.spoofTimesToOneSession(trials,{'imageTime',timefieldname});
            imageTimes = eventTimes.imageTime;
            eventTimes = eventTimes.(timefieldname);
            
            if(~isempty(contingency))
                contingentEventTimes = ImageHelper.spoofTimesToOneSession(trials,contingency(1));
                contingentEventTimes = contingentEventTimes.(contingency{1});
                contingentEventTimeRange = contingency{2};
                contingentEventKeep = contingency{3};
            end
            
            for i=1:length(eventTimes)
                if(~isempty(contingency))
                    % Ensure this event time matches the given contingency,
                    % otherwise skip it (and ultimately delete this image)
                    contingencyDeltaT = contingentEventTimes - eventTimes(i);
                    meetContingency = any(contingencyDeltaT > contingentEventTimeRange(1) & contingencyDeltaT < contingentEventTimeRange(2));
                    if(contingentEventKeep == meetContingency)
                        imdelete(1,1,1,i) = true;
                        %disp(sprintf('Event #%d does not match contingency.',i));
                        continue;
                        
                    end
                end
                
                
                % find nearest image
                [closestIx] = find(imageTimes > eventTimes(i),1,'first');
                
                % Go backwards
                currentImageIndex = closestIx-1;
                while(currentImageIndex > 0                     && (eventTimes(i) - imageTimes(currentImageIndex)) < parseLength)
                    frameOffset = floor((imageTimes(currentImageIndex) - eventTimes(i))/FRAMETIME) + ix0;
                    if(frameOffset > 0 && frameOffset < numImages)
                        im(:,:,frameOffset,i) = images(:,:,currentImageIndex);
                    end
                    currentImageIndex = currentImageIndex-1;
                end
                
                % Go forwards
                currentImageIndex = closestIx;
                while(currentImageIndex <= length(imageTimes)   && (imageTimes(currentImageIndex) - eventTimes(i)) < parseLength)
                    frameOffset = floor((imageTimes(currentImageIndex) - eventTimes(i))/FRAMETIME) + ix0;
                    if(frameOffset > 0 && frameOffset < numImages)
                        im(:,:,frameOffset,i) = images(:,:,currentImageIndex);
                    end
                    currentImageIndex = currentImageIndex+1;
                end
            end
            
            im(:,:,:,imdelete) = [];
        end
        
        function [pre,ix0,post] = alignImages(trials,images,timefieldname)
            alignmentTimes = {trials.(timefieldname)};
            pre = cell(length(trials),1);
            post = cell(length(trials),1);
            ix0 = 1;
            
            for i=1:length(trials)
                [pre{i},post{i}] = ImageHelper.getPreAndPostImages(images(:,:,ImageHelper.imageIndices(trials(i))),trials(i).imageTime,alignmentTimes{i});
                ix0 = max(ix0,size(pre{i},3)+1);
            end
            
            
            pre = ImageHelper.catPreImages(pre);
            post = ImageHelper.catPostImages(post);
            
            if(nargout == 2)
                pre = ImageHelper.catPreToPost(pre,post);
            end
        end
       
        function [preimages,postimages] = getPreAndPostImages(images,imageTimes,alignmentTime)
            badData = false;
            if(size(images,3) ~= numel(imageTimes))
                warning('Size mismatch between images (%d) and imageTimes (%d) - returning empty datasets.',size(images,3),numel(imageTimes));
                badData = true;
            end
            if(any(diff(imageTimes) < 0))
                warning('imageTimes must be monotonically increasing - returning empty datasets.');
                badData = true;
            end
            if(isempty(alignmentTime))
                %warning('alignmentTime is empty - returning empty datasets.');
                badData = true;
            end
            if(badData)
                preimages = single(zeros(size(images,1),size(images,2), 0));
                postimages = single(zeros(size(images,1),size(images,2), 0));
                return;
            end
            
            splitIx = find(imageTimes > alignmentTime,1,'first');
            
            if(isempty(splitIx))
                % no image time is greater than alignment time = all images
                % occurred in preimage period
                preimages = images;
                postimages = single(zeros(size(images,1),size(images,2), 0));
            elseif(splitIx == 1)
                % all image times are greater than alignment time = all
                % images occurred in postimage period
                preimages = single(zeros(size(images,1),size(images,2), 0));
                postimages = images;
            else
                % splitIx is somewhere between 2 and the last image
                preimages = images(:,:,1:(splitIx-1));
                postimages = images(:,:,splitIx:end);
                
            end
        end
        
        function I = catPreToPost(preimages,postimages)
            I = cat(3,preimages,postimages);
        end
        
        function I = catPreImages(imcell)
            %{
            Preimages are images before the alignment point.
            They are concatenated right-aligned with left-padding.
            %}
            
            SM = [0 0 0];
            for i=1:length(imcell)
                %S = sizen(imcell{i},3);
                S = size(imcell{i});
                SM = max(SM,S);
            end
            
            I = nan([SM(1) SM(2) SM(3) length(imcell)],'single');
            for i=1:length(imcell)
                %S = sizen(imcell{i},3);
                S = size(imcell{i});
                I(:,:,end-S(3)+1:end,i) = imcell{i};
                %I(:,:,end-S+1:end,i) = imcell{i};
            end
        end
        
        function I = catPostImages(imcell)
            %{
            Postimages are images before the alignment point.
            They are concatenated left-aligned with right-padding.
            %}
            SM = [0 0 0];
            for i=1:length(imcell)
                %S = sizen(imcell{i},3);
                S = size(imcell{i});
                SM = max(SM,S);
            end
            
            I = nan([SM(1) SM(2) SM(3) length(imcell)],'single');
            for i=1:length(imcell)
                %S = sizen(imcell{i},3);
                S = size(imcell{i});
                I(:,:,1:S(3),i) = imcell{i};
            end
        end
            
        %% Category methods
        % Sometimes you want to consider subsets of trials. These methods
        % use user defined "categories" (1 to N) to manipulate and sum
        % subsets of trials
        function [count, sum1, sum2, themean, thestd, thesterr] = sumImagesByCategory(alignedimages,category)
%             if(ndims(alignedimages) ~= 4)
%                 error('4D data expected for alignedimages.');
%             end
            
            numcat = max(category);
            [s1,s2,s3,s4] = size(alignedimages);
            
            count = zeros(s1,s2,s3,numcat);
            sum1  = zeros(s1,s2,s3,numcat);
            sum2  = zeros(s1,s2,s3,numcat);
            imageDataHolder = nan(s1,s2,s3,numcat);
            
%             for i=1:numcat
%                 fprintf('Counting category %d/%d... \n',i,numcat);
%                 im = alignedimages(:,:,:,category == i);
%                 count(:,:,:,i) = sum(~isnan(im),4);
%                 sum1(:,:,:,i) = nansum(im,4);
%                 if(nargout >= 3)
%                     sum2(:,:,:,i) = nansum(im.^2,4);
%                 end
%             end
            
            for i=1:size(alignedimages,4)
%                 imageDataHolder(:) = nan;
%                 imageDataHolder(:,:,:,category(i)) = alignedimages(:,:,:,i);
%                 count(~isnan(imageDataHolder)) = count(~isnan(imageDataHolder))+1;
%                 sum1 = nansum(cat(5,sum1,imageDataHolder),5);
%                 sum2 = nansum(cat(5,sum2,imageDataHolder.^2),5);
                
                valid_pixels = find(~isnan(alignedimages(:,:,:,i)));
                [d1,d2,d3]   = ind2sub([s1 s2 s3],valid_pixels);
                
                write_from   = sub2ind([s1 s2 s3 s4],d1,d2,d3,i*ones(size(d1)));
                write_to     = sub2ind([s1 s2 s3 numcat],d1,d2,d3,category(i)*ones(size(d1)));
                
                count(write_to) = count(write_to)+1;
                sum1(write_to)  = sum1(write_to) + alignedimages(write_from);
                sum2(write_to)  = sum2(write_to) + alignedimages(write_from).^2;
            end
            
            themean = sum1./count;
            thestd = sqrt( (sum2 - (sum1.^2)./count) ./ (count-1) );
            thesterr = thestd ./ sqrt(count);
        end
        function [Icell,Tcell] = seperateCategorySubsets(I,T,category)
            numcategories = max(category);
            Tcell = cell(numcategories,1);
            Icell = cell(numcategories,1);
            for cat=1:numcategories
                Tcell{cat} = T(category == cat);
                
                % Count images
                imcount = 0;
                for tri = 1:length(Tcell{cat})
                    imcount = imcount + length(Tcell{cat}(tri).imageTime);
                end
                
                II = zeros(size(I,1),size(I,2),imcount);
                imcount = 0;
                
                for tri = 1:length(Tcell{cat})
                    imix = ImageHelper.imageIndices(Tcell{cat}(tri));
                    II(:,:,imcount+1:imcount+length(imix)) = I(:,:,imix);
                    Tcell{cat}(tri).imageData = imcount+1;
                    imcount = imcount+length(imix);
                end
                Icell{cat} = II;
            end
        end
        function [I,T] = combineCategorySubsets(Icell,Tcell)
            % Count and pre-allocate T and I
            numIm = 0;
            numTr = 0;
            for cat=1:length(Tcell)
                numTr = numTr + length(Tcell{cat});
                numIm = numIm + size(Icell{cat},3);
            end
            [s1,s2,~,~] = size(Icell{1});
            
            % Preallocate
            I = zeros(s1,s2,numIm,'single');
            T(numTr) = Tcell{1}(1);
            tcount = 0;
            icount = 0;
            
            for cat=1:length(Tcell)
                numTr = length(Tcell{cat});
                % Update imageData indicies
                for tri=1:numTr
                    Tcell{cat}(tri).imageData = Tcell{cat}(tri).imageData + icount;
                end
                T(tcount+1:tcount+numTr) = Tcell{cat};
                tcount = tcount+numTr;
                
                numIm = size(Icell{cat},3);
                I(:,:,icount+1:icount+numIm) = Icell{cat};
                icount = icount+numIm;
            end
        end
        
        function cat = categorize(T,catstr)
            % Categories is a string which indexes a field of T
            catData = cell(length(T),1);
            for i=1:length(T)
                eval(['catData{i} = T(i).' catstr ';']);
            end
            [~,~,cat] = unique(catData,'stable');
        end
        function catI = categorizeTrialToImage(T,catT)
            % categorize works on the basis of trials. This converts trial
            % categories into image categories.
            numImg = T(end).imageData + length(T(end).imageTime) - 1;
            
            catI = nan(numImg,1);
            for i=1:length(T)
                catI(ImageHelper.imageIndices(T(i))) = catT(i);
            end
        end
        
        function [count, sum1, sum2, themean, thestd, thesterr] = combineCategorySums(count, sum1, sum2,d1,d2,d3,d4)
            if(~isempty(d1))
                count = nansum(count(d1,:,:,:),1);
                sum1 = nansum(sum1(d1,:,:,:),1);
                sum2 = nansum(sum2(d1,:,:,:),1);
            end
            if(~isempty(d2))
                count = nansum(count(:,d2,:,:),2);
                sum1 = nansum(sum1(:,d2,:,:),2);
                sum2 = nansum(sum2(:,d2,:,:),2);
            end
            if(~isempty(d3))
                count = nansum(count(:,:,d3,:),3);
                sum1 = nansum(sum1(:,:,d3,:),3);
                sum2 = nansum(sum2(:,:,d3,:),3);
            end
            if(~isempty(d4))
                count = nansum(count(:,:,:,d4),4);
                sum1 = nansum(sum1(:,:,:,d4),4);
                sum2 = nansum(sum2(:,:,:,d4),4);
            end
            
            themean = sum1./count;
            thestd = sqrt( (sum2 - (sum1.^2)./count) ./ (count-1) );
            thesterr = thestd ./ sqrt(count);
        end
        %% Nusiance Regression
        
        function I2 = downsampleImages(I)
            % Always works on factors of 2
            [D1,D2,D3,D4] = size(I);            
            I2 = (I(1:2:end,1:2:end,:,:) + I(2:2:end,1:2:end,:,:) + I(1:2:end,2:2:end,:,:) + I(2:2:end,2:2:end,:,:))/4;
        end
        function I2 = upsampleImages(I)
            % Always works on factors of 2
            [D1,D2,D3,D4] = size(I);
            I2 = zeros(2*D1,2*D2,D3,D4);
            I2(1:2:end,1:2:end,:,:) = I;
            I2(2:2:end,1:2:end,:,:) = I;
            I2(1:2:end,2:2:end,:,:) = I;
            I2(2:2:end,2:2:end,:,:) = I;
        end
        
        function I2 = globalMeanTimecourse(I)
            I2 = mean(mean(I,2),1);
        end
        
        % Regress out polynomials or other functions
        function regmat = getPolynomialRegressor(N,order)
            if(nargin == 0)
                order = 2;
            end
            
            regmat = zeros(N,order);
            for i=0:order
                regressor = linspace(0,1,N).^i;
                regmat(:,i+1) = regressor';
            end
        end
        function [I,b] = regressFromImages(I,varargin)
            nvarargs = length(varargin);
            [X,Y,N,C] = size(I);
            
            const_regmat = zeros(size(I,3),0);
            for i=1:nvarargs
                if(ismatrix(varargin{i})) 
                    % 2-D is assumed to be a component of the regression matrix that is constant for all pixels
                    const_regmat = cat(2,const_regmat,varargin{i});
                end
            end
            
            b = nan(X,Y,1);
            for x=1:X
                for y=1:Y
                    % Add per-pixel regressors
                    regmat = const_regmat;
                    for i=1:nvarargs
                        if(~ismatrix(varargin{i}))
                            regmat = cat(2,regmat,squeeze(varargin{i}(x,y,:)));
                        end
                    end
                    
                    for c=1:C
                        ts = squeeze(I(x,y,:,c));
                        if(~all(isnan(ts)))
                            b(x,y,1:size(regmat,2)) = regress(ts,regmat);
                            I(x,y,:,c) = permute(ts - regmat*reshape(b(x,y,:),[],1),[2 3 1 4]);
                        end
                    end
                end
            end
        end
        
        %% Eye Data Methods
        % If you want the raw eye data for saccade detection or other
        % high-frequency analysis, use this.
        function [X,Y] = convertEyeData(T)
            Xraw = double(T.eyeData(1:2:end));
            Yraw = double(T.eyeData(2:2:end));
            
            X = T.eyeDataCalibration.cal.m11 * Xraw   +   T.eyeDataCalibration.cal.m21 * Yraw   +   T.eyeDataCalibration.cal.tX;
            Y = T.eyeDataCalibration.cal.m12 * Xraw   +   T.eyeDataCalibration.cal.m22 * Yraw   +   T.eyeDataCalibration.cal.tY;
        end
        
        
        function [breakTime] = checkForFixationBreak(T,window)
            [X,Y] = ImageHelper.convertEyeData(T);
            eyeTime = linspace(double(T.trialStartTime),double(T.trialEndTime(1)),length(X));
            
            fixTime  = double(T.prestimulusTime);
            stimTime = double(T.stimulusOnTime);
            behTime = double(ImageHelper.getEndOfFixationPeriod(T));
            
            baseline = find(eyeTime > fixTime & eyeTime < stimTime,50,'last');
            baseX = mean(X(baseline));
            baseY = mean(Y(baseline));
            
            X(eyeTime < stimTime) = [];
            Y(eyeTime < stimTime) = [];
            eyeTime(eyeTime < stimTime) = [];

            bt =     (abs(X-baseX) > window | ...
                      abs(Y-baseY) > window );
            bt = find((bt + circshift(bt,[0 -1]) + circshift(bt,[0 -2]) + ...
                            circshift(bt,[0 -3]) + circshift(bt,[0 -4]) + ...
                            circshift(bt,[0 -5]) + circshift(bt,[0 -6]) + ...
                            circshift(bt,[0 -7]) + circshift(bt,[0 -8]) +...
                            circshift(bt,[0 -9]) + circshift(bt,[0 -10]) ...
                            ) == 11,1,'first');
            
            breakTime = eyeTime(bt);
            if(isempty(breakTime))
                breakTime = Inf;
            end;
        end
        
        
        % Otherwise, use this to compare with images. NOW CIRCSHIFTED!
        function [eyeTime,X,Y,Xp,Yp] = convertEyeDataToImageTimes(T,eyeImageTimeOffset,eyeSamplesPerImage)
            eyeSampleTime = ImageHelper.defaultEyeSampleTime;
            if(nargin == 1)
                eyeImageTimeOffset = 0;
                eyeSamplesPerImage = 1;
            end
            
            % Real units
            [X,Y] = ImageHelper.convertEyeData(T);
            
            % Correct number of samples, fill in any gaps, assign a time to
            % each sample
            [X,Y,eyeTime] = ImageHelper.fillTimePeriodWithEyeData(X,Y, ...
                double(T.trialStartTime)/100,double(T.trialEndTime)/100,eyeSampleTime);
            
            % Create a struture to "match" the images
            [X,Y,eyeTime] = ImageHelper.groupAndAverageEyeDataByTime(X,Y,eyeTime,double(T.imageTime)/100,eyeImageTimeOffset,eyeSamplesPerImage);
        end
        
        % These functions return a meaningful eye data structure for future
        % analyses
        function [X,Y,ix] = alignEyeData(T,timefieldname)
            % Returns a TxN matrix of eye-positions per trial, aligned to
            % some timepoint.
            eyeSampleTime = 1000/240;
            
            % As in other functions, we assume that eye data span the
            % interval from TrialStart to TrialEnd with no interruptions.
            % However, we only accept this as valid if there is not more
            % than 300 ms of "missing" eye data.
            Xs_pre = cell(size(T));
            Ys_pre = cell(size(T));
            Xs_pos = cell(size(T));
            Ys_pos = cell(size(T));
            
            for i=1:length(T)
                missingDataTime = ImageHelper.timeOfMissingEyeData(T(i));
                alignTime = T(i).(timefieldname)/100;
                
                if(missingDataTime < 300 && ~isempty(alignTime))
                    [X,Y] = ImageHelper.convertEyeData(T(i));
                    [X,Y,eyeDataTime] = ImageHelper.fillTimePeriodWithEyeData(X,Y,double(T(i).trialStartTime)/100,double(T(i).trialEndTime)/100,eyeSampleTime);

                    transitionIndex = find(eyeDataTime > alignTime,1,'first');
                    if(~isempty(transitionIndex))
                        if(transitionIndex > 1)
                            Xs_pre{i} = X(1:(transitionIndex-1));
                            Ys_pre{i} = Y(1:(transitionIndex-1));
                        end
                        Xs_pos{i} = X(transitionIndex:end);
                        Ys_pos{i} = Y(transitionIndex:end);
                    end
                end
            end
            
            longest_pre= 0;
            longest_pos= 0;
            for i=1:length(T)
                longest_pre = max(length(Xs_pre{i}),longest_pre);
                longest_pos = max(length(Xs_pos{i}),longest_pos);
            end            
            
            ix = longest_pre+1;
            X = nan(length(T),longest_pre+longest_pos);
            Y = nan(length(T),longest_pre+longest_pos);
            
            for i=1:length(T)
                X(i,(ix-length(Xs_pre{i})):(ix-1)) = Xs_pre{i};
                Y(i,(ix-length(Ys_pre{i})):(ix-1)) = Ys_pre{i};
                X(i,ix-1+(1:length(Xs_pos{i}))) = Xs_pos{i};
                Y(i,ix-1+(1:length(Ys_pos{i}))) = Ys_pos{i};
            end
        end
        
        % The Big Kahuna - This returns a HEAVILY MODIFIED I and T dataset
        % Only images with viable eye-signals are retained.
        function [I2,T,regmat,beta] = regressOutEyeMovements(I,T,realRegression)
            if(nargin==2)
                realRegression = true;
            end
            
            convertToSingle = strcmp(class(I),'single');
            
            disp('Pre-allocating.');
            I2 = zeros(size(I));
            bigX = zeros(size(I,3),1);
            bigY = zeros(size(I,3),1);
            imcounter = 0;
            
            x_counter = 0;
            y_counter = 0;
            i_counter = 0;
            
            disp('Gathering images and eye positions.');
            for i=1:length(T)
                if(isempty(T(i).eyeData))
                    T(i).imageData = imcounter+1;
                    T(i).imageTime = [];
                    continue;
                end
                
                % First grab raw eyeData from just after stimulus onset
                stimOnEyeDataIndex = floor(double((T(i).stimulusOnTime - T(i).trialStartTime))/100/4.1667);
                % Start on X
                if(mod(stimOnEyeDataIndex,2)==1)
                    stimOnEyeDataIndex = stimOnEyeDataIndex + 1;
                end
                eyeRange1 = stimOnEyeDataIndex + 120;
                eyeRange2 = stimOnEyeDataIndex + 240;
                
                [X,Y] = ImageHelper.convertEyeData(T(i));
                
                if(eyeRange1 < length(T(i).eyeData))
                    xx = X( eyeRange1:min(eyeRange2,length(X)));
                    yy = Y( eyeRange1:min(eyeRange2,length(Y)));
                end
               
                x_counter = x_counter + sum(xx);
                y_counter = y_counter + sum(yy);
                i_counter = i_counter + length(xx);
                
                [~,X,Y] = ImageHelper.convertEyeDataToImageTimes(T(i));
                II = double(I(:,:,ImageHelper.imageIndices(T(i))));
                
                if(isempty(II))
                    T(i).imageData = imcounter+1;
                    T(i).imageTime = [];
                    continue;
                end
                
                X = circshift(X,[2 0]);
                Y = circshift(Y,[2 0]);
                
                delix = isnan(X) | isnan(Y);
                X(delix) = [];
                Y(delix) = [];
                II(:,:,delix) = [];
                
                % Demean eye data to prevent differences in eye calibration
                % error from having undue influence across trials
                X = demean(X);
                Y = demean(Y);
                
                T(i).imageData = imcounter+1;
                T(i).imageTime(delix) = [];
                
                len = length(X);
                bigX(imcounter+1:imcounter+len) = X;
                bigY(imcounter+1:imcounter+len) = Y;
                I2(:,:,imcounter+1:imcounter+len) = II;
                imcounter = imcounter+len;
            end
            
            I2 = I2(:,:,1:imcounter);
            bigX = bigX(1:imcounter);
            bigY = bigY(1:imcounter);
            big1 = ones(size(bigX));
            
            if(realRegression)
                regmat = [bigX,bigY,big1];
            else
                regmat = big1;
            end
            
            fprintf('Mean Stim-On Eye Position: (%f, %f)\n',x_counter/i_counter,y_counter/i_counter);
            
            [I2,beta] = ImageHelper.regressFromImages(I2,regmat);
            if(convertToSingle)
                I2 = single(I2);
            end
        end
        function [I,T,R] = regressOutEyeMovementsByCategory(I,T,C,realRegression)
            %{
                Helper function to sort images by category before
                regression. Useful when categories represent distinct
                visual stimuli.
            %}
            if(nargin == 3)
                realRegression = true;
            end
            
            [Icell,Tcell] = ImageHelper.seperateCategorySubsets(I,T,C);
            R = cell(9,1);
            
            
            for i=1:length(Tcell)
                [Icell{i},Tcell{i},R{i}] = ImageHelper.regressOutEyeMovements(Icell{i},Tcell{i},realRegression);
            end
            
            [I,T] = ImageHelper.combineCategorySubsets(Icell,Tcell);
        end
        
        % Remainder of eye data methods generally serve the above method.
        function [X,Y] = downsampleEyeData(Xin,Yin,factor)
            len = length(Xin);
            X = zeros(ceil(len/factor),1);
            Y = zeros(ceil(len/factor),1);
            i1 = 1;
            i2 = 1;
            while i1 <= len
                if(i1+factor-1 > len)
                    X(i2) = mean(Xin(i1:end));
                    Y(i2) = mean(Yin(i1:end));
                else
                    X(i2) = mean(Xin(i1:i1+factor-1));
                    Y(i2) = mean(Yin(i1:i1+factor-1));
                end
                
                i1 = i1+factor;
                i2 = i2+1;
                
            end
        end
        function [X,Y] = interpolateEyeData(Xin,Yin,targetLength)
            lenE = length(Xin);
            X = interp1(linspace(0,1,lenE),Xin,linspace(0,1,targetLength));
            Y = interp1(linspace(0,1,lenE),Yin,linspace(0,1,targetLength));
        end
        function [X,Y,eyeDataTimes] = fillTimePeriodWithEyeData(Xin,Yin,timeOn,timeOff,eyeSamplePeriod)
            if(isempty(timeOn) || isempty(timeOff) || length(timeOn) > 1 || length(timeOff) > 1)
                X = [];
                Y = [];
                eyeDataTimes = [];
                return;
            end

            targetLength = ceil((timeOff-timeOn)/eyeSamplePeriod);
            [X,Y] = ImageHelper.interpolateEyeData(Xin,Yin,targetLength);
            eyeDataTimes = linspace(timeOn,timeOff,targetLength);
        end
        function [X,Y,eyeTime] = groupAndAverageEyeDataByTime(Xin,Yin,eyeTime,averageTime,eyeImageTimeOffset,eyeSamplesPerImage)
            binWidth = 200; % ms
            
            bin_end = averageTime;
            bin_begin = averageTime - binWidth;
            
            lenout = length(bin_end);
            
            X = nan(lenout,1);
            Y = nan(lenout,1);
            
            for i=1:lenout
                inTimeBin = eyeTime >= bin_begin(i) & eyeTime <= bin_end(i);
                X(i) = mean(Xin(inTimeBin));
                Y(i) = mean(Yin(inTimeBin));
            end
        end      
        function missingTime = timeOfMissingEyeData(T,eyeSampleTime)
            if(nargin==1)
                eyeSampleTime = ImageHelper.defaultEyeSampleTime;
            end
            startTime = double(T.trialStartTime);
            endTime   = double(T.trialEndTime);
            estimatedSamples = (endTime-startTime)/100/eyeSampleTime;
            actualSamples = length(T.eyeData)/2;
            missingTime = abs(estimatedSamples - actualSamples)*eyeSampleTime;
        end
        function t = defaultEyeSampleTime()
            t = 1000/240;
        end
        %% Full one-stop nusiance regression
        function I2 = nuisanceRegressionOneTrial(I,T,eyeSampleRate,cameraSampleRate,globalSignalMask,polynomialOrder)
            if(nargin == 4 || isempty(globalSignalMask))
                globalSignal = ImageHelper.globalMeanTimecourse(I);
            else
                if(isvector(globalSignalMask))
                    globalSignal = mean(mean(I(globalSignalMask,globalSignalMask,:),2),1);
                else
                    globalSignal = mean(mean(I(repmat(globalSignalMask,[1 1 size(I,3)])),2),1);
                end
            end
            
            imageCount = size(I,3);
            [trialEyeX,trialEyeY] = ImageHelper.convertEyeData(T);
            [trialEyeX,trialEyeY] = ImageHelper.interpolateEyeData(trialEyeX,trialEyeY,imageCount/cameraSampleRate*eyeSampleRate);
            [trialEyeX,trialEyeY] = ImageHelper.downsampleEyeData(trialEyeX,trialEyeY,eyeSampleRate/cameraSampleRate);
            
            poly = ImageHelper.getPolynomialRegressor(imageCount,polynomialOrder);
            
            I2 = ImageHelper.regressFromImages(I,[poly squeeze(globalSignal) circshift([trialEyeX trialEyeY],[2 0])]);
        end
        function I = nuisanceRegression(I,T,F,globalSignalMask)
            if(nargin == 3)
                globalSignalMask = [];
            end
            
            imageDataDeviceIndex = find(strcmp('imageData',F.dataGroupParam.singleData.dataName),1,'first');
            eyeDataDeviceIndex = find(strcmp('eyeData',F.dataGroupParam.singleData.dataName),1,'first');
            imageRate = round(1000/F.dataGroupParam.singleData.timing(imageDataDeviceIndex));
            eyeRate = round(1000/F.dataGroupParam.singleData.timing(eyeDataDeviceIndex));
           
            fprintf(1,'Regressing nuisance variables... Trials Done: 000000/000000');
            numTrials = length(T);
            for i=1:numTrials
                imgix = T(i).imageData : T(i).imageData+length(T(i).imageTime)-1;
                I(:,:,imgix) = ImageHelper.nuisanceRegressionOneTrial( ...
                    I(:,:,imgix),T(i),eyeRate,imageRate,globalSignalMask,0);
                fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b%06d/%06d',i,numTrials);
            end
            
        end
        
        %% Make it easier to index into trials
        function date = getDateModified(T)
            date = cell(length(T),1);
            for i=1:length(T)
                info = dir(T(i).fullDataPath);
                longdate = info.date;
                date{i} = longdate(1:find(longdate==' ',1,'first')-1);
            end
        end
        
        function imix = imageIndices(T)
            imix =  T.imageData : T.imageData+length(T.imageTime)-1;
        end
        function imix = firstPostStimImageIndex(T)
            
            imix = find(T.imageTime > T.stimulusOnTime,1,'first');
            ix = ImageHelper.imageIndices(T);
            imix = ix(imix);
        end
        
        function [I,T] = retainPrebehaviorImages(I,T,graceTimeMS,breakWindow)
            if(nargin == 3)
                breakWindow = inf;
            end
            imcounter = 0;
            imremove = false(size(I,3),1);
            grace = graceTimeMS*100;
            
            for tri = 1:length(T)
                behaviorTime = min(ImageHelper.getEndOfFixationPeriod(T(tri),breakWindow));
                if(isempty(behaviorTime))
                    % Incomplete trial
                    behaviorTime = -grace;
                end
                
                badImages = find(T(tri).imageTime > (behaviorTime+grace));
                
                imix = ImageHelper.imageIndices(T(tri));
                imremove(imix(badImages)) = true;
                
                T(tri).imageTime(badImages) = [];
                
                T(tri).imageData = imcounter+1;
                imcounter = imcounter+length(T(tri).imageTime);
            end
            
            I(:,:,imremove) = [];
        end
        function [I,T] = retainPostbehaviorImages(I,T,graceTimeMS,breakWindow)
            if(nargin == 3)
                breakWindow = inf;
            end
            imcounter = 0;
            imremove = false(size(I,3),1);
            grace = graceTimeMS*100;
            
            for tri = 1:length(T)
                behaviorTime = ImageHelper.getEndOfFixationPeriod(T(tri),breakWindow);
                if(isempty(behaviorTime))
                    % Incomplete trial
                    behaviorTime = -grace;
                end
                
                badImages = find(T(tri).imageTime < (behaviorTime-grace));
                
                imix = ImageHelper.imageIndices(T(tri));
                imremove(imix(badImages)) = true;
                
                T(tri).imageTime(badImages) = [];
                
                T(tri).imageData = imcounter+1;
                imcounter = imcounter+length(T(tri).imageTime);
            end
            
            I(:,:,imremove) = [];
        end
        function [I,T] = removePrePrestimulusImages(I,T,graceTimeMS)
            imcounter = 0;
            imremove = false(size(I,3),1);
            grace = graceTimeMS*100;
            
            for tri = 1:length(T)
                cutoffTime = T(tri).prestimulusTime;
                if(isempty(cutoffTime))
                    % Incomplete trial
                    cutoffTime = Inf;
                end
                badImages = find(T(tri).imageTime < (cutoffTime-grace));
                
                imix = ImageHelper.imageIndices(T(tri));
                imremove(imix(badImages)) = true;
                
                T(tri).imageTime(badImages) = [];
                
                T(tri).imageData = imcounter+1;
                imcounter = imcounter+length(T(tri).imageTime);
                
            end
            
            I(:,:,imremove) = [];
        end
        
        function [I,T] = assertOverTrials(I,T,assertion)
            % 1 - Find trials to keep
            trialkeep = false(size(T));
            imagekeep = false(size(I,3),1);
            imcount = 0;
            for t=1:length(T)
                if assertion(T(t))
                    trialkeep(t) = true;
                    trialImageLength = length(T(t).imageTime);
                    imagekeep(imcount+1:imcount+trialImageLength) = true;
                    T(t).imageData = imcount+1;
                    imcount = imcount+trialImageLength;
                end
            end
            
            T(~trialkeep) = [];
            I(:,:,~imagekeep) = [];
        end
        
        function fixateOffTime = getEndOfFixationPeriod(trial,breakWindow)
            if(nargin==2)
                fixateOffTime = min([trial.trialEndTime trial.deviantOnTime trial.saccadeTime trial.stimulusOffTime ImageHelper.checkForFixationBreak(trial,breakWindow)]);
            else
                fixateOffTime = min([trial.trialEndTime trial.deviantOnTime trial.saccadeTime trial.stimulusOffTime]);
            end
        end
        
        %% Unit transformations
        function [I,T] = convertToPercentChangeFromBaseline(I,T)
            % NOTE - Trials without at least one post-stimulus image
            % (trials which return nothing for firstPostStimImageIndex)
            % have *all* images removed, as their units in terms of 
            % percent change are undefined.
            %
            % Image values are set to NaN under the following conditions
            % 1 - firstPostStimImageIndex returns []
            % 2 - in a given pixel, ANY value in the trial reaches 0 or
            %     65535
            % 3 - There must be at least 1 image PRIOR to
            %     firstPostStimImageIndex

            [X,Y,L] = size(I);
            
            for t = 1:length(T)                
                imix = ImageHelper.imageIndices(T(t));
                divIx = ImageHelper.firstPostStimImageIndex(T(t));
                
                if(isempty(imix) || isempty(divIx) || divIx == imix(1))
                    % Criterion 1 or 3 reached
                    I(:,:,imix) = nan;
                    continue;
                end
                
                images = I(:,:,imix);
                
                void = any(images == 0 ,3);
                
                images(repmat(void,[1 1 size(images,3)])) = nan;
                
                %I(:,:,imix) = images ./ repmat(nanmean(I(:,:,imix(1):divIx),3),[1 1 length(imix)]);
                I(:,:,imix) = images ./ repmat(nanmean(I(:,:,divIx),3),[1 1 length(imix)]);
                
                %I(:,:,imix) = images ./ repmat(nanmean(I(:,:,divIx+(-4:0)),3),[1 1 length(imix)]);
                
                %I(:,:,imix) = images;
                %I(:,:,imix) = images ./ repmat(mean(I(:,:,imix),3),[1 1 length(imix)]);
            end
        end
        
        %% Image Affine Transformations
        function IA = affineTransform(I,A,border)
            if(nargin==2)
                border = 0;
            end
            
            disp('Pre-allocating coordinate system and output structure (WATCH TOP!).');
            IA = zeros(size(I,1)+2*border, size(I,2)+2*border, size(I,3),'uint16');
            
            
            [X2,Y2] = meshgrid(1:size(I,2)+border*2, 1:size(I,1)+border*2   );
            X1 = X2(border+1:end-border,border+1:end-border);
            Y1  =Y2(border+1:end-border,border+1:end-border);
            
            I2 = ones(size(X2));
            
            targetPixelCoordinates = [X2(:)'; Y2(:)'; I2(:)'];
            target = IA(:,:,1);
            
            for i=1:size(I,3)
                if(isnan(A(1,1,i)))
                    IA(1:5,:,i) = ImageHelper.BAD_IMAGE_INT5;
                    continue;
                end
                
                % Syntax reminder: this is inv(A)*targetPixelCoordinates
                R = A(:,:,i)\targetPixelCoordinates;
                
                X = reshape(R(1,:),size(target));
                Y = reshape(R(2,:),size(target));
                
                if(isa(I,'single'))
                    IA(:,:,i) = uint16(interp2(X1,Y1,I(:,:,i),X,Y)        );
                else
                    IA(:,:,i) = uint16(interp2(X1,Y1,single(I(:,:,i)),X,Y));
                end
            end
        end
        
        function bad = badImages(I)
            bad = false(1,1,size(I,3));
            for i=1:size(I,3)
                if(all(I(1:5,:,i) == ImageHelper.BAD_IMAGE_INT5))
                    bad(i) = true;
                end
            end
        end
        
        function A = fillInBadTransforms(A)
            % A bad transformation is one with all NANs. Program will find
            % the two flanking good transforms and average them.
            
            for i=1:size(A,3)
                if(isnan(A(:,:,i)))
                    down = i;
                    while(isnan(A(:,:,down)))
                        down = down-1;
                        if(down < 1)
                            break;
                        end
                    end
                    up = i;
                    while(isnan(A(:,:,up)))
                        up = up+1;
                        if(up > size(A,3))
                            break;
                        end
                    end

                    if(up > size(A,3) && down < 1)
                        error('No good transformations found!');
                    elseif(up > size(A,3))
                        A(:,:,i) = A(:,:,down);
                    elseif(down < 1)
                        A(:,:,i) = A(:,:,up);
                    else
                        A(:,:,i) = A(:,:,up)/2+A(:,:,down)/2;
                    end
                end
            end
            
        end
        
        function I = reorientImages(I)
            % Convert images to be oriented upright with real L/R on image
            % L/R
            for i=1:size(I,3)
                for j=1:size(I,4)
                    for k=1:size(I,5)
                        I(:,:,i,j,k) = fliplr(I(:,:,i,j,k)');
                    end
                end
            end
        end
        
        %% Statistics helpers
        function activation = getActivationForTrials(T,I,low,high)
            % Returns the average difference image between low and high
            % indices relative to stimulus onset. Any images not present
            % are filled in with NAN.
            activation = nan(size(I,1),size(I,2),length(T));
            for tri = 1:length(T)
                imix = ImageHelper.imageIndices(T(tri));
                stimOn = ImageHelper.firstPostStimImageIndex(T(tri));
                
                if(isempty(stimOn))
                    continue;
                end
                
                thisLow  = stimOn+low;
                thisHigh = stimOn+high;
                
                if(~all(ismember(thisLow,imix)) || ~all(ismember(thisHigh,imix)))
                    continue;
                end
                
                activation(:,:,tri) = nanmean(I(:,:,thisHigh),3)-nanmean(I(:,:,thisLow),3);
            end
        end
        
        function spoofedTimes = spoofTimesToOneSession(T,fieldNames,fileOffsetTime)
            % times = spoofTimesToOneSession(T,fieldNames)
            % T = Trial Struct
            % fieldNames = cell array of strings
            % times = cell array of times
            
            if(nargin==2)
                fileOffsetTime = 10000000; % 100 seconds
            end
            if(ischar(fieldNames))
                fn = fieldNames;
                fieldNames = cell(1,1);
                fieldNames{1}  =fn;
            end
            
            % T-zero = first trial start time
            firstTrialTimeReal = double(T(1).trialStartTime);
            firstTrialTimeFake = 0;
            lastTrialStartTime = 0;
            currentDataPath = T(1).fullDataPath;
            
            for j=1:length(fieldNames)
                command = sprintf('spoofedTimes.%s = [];',fieldNames{j});
                eval(command);
            end
            for i=1:length(T)
                if(~strcmp(T(i).fullDataPath,currentDataPath))
                    firstTrialTimeReal = double(T(i).trialStartTime);
                    firstTrialTimeFake = lastTrialStartTime + fileOffsetTime;
                    currentDataPath = T(i).fullDataPath;
                end
                
                for j=1:length(fieldNames)
                    command = (sprintf('spoofedTimes.%1$s(end+1:end+length(T(i).%1$s)) = double(T(i).%1$s) - firstTrialTimeReal + firstTrialTimeFake;',fieldNames{j}));
                    eval(command);
                end
                lastTrialStartTime = double(T(i).trialStartTime) - firstTrialTimeReal + firstTrialTimeFake;
            end
            
        end
        
        %% Debugging
        function bool = checkTrialImageIntegrity(T)
            if(~iscell(T))
                Tcell{1} = T;
            else
                Tcell = T;
            end
            
            bool = true;
            for j=1:length(Tcell)
                T = Tcell{j};
                
                largestImageIndex = -1;
                for i=1:length(T)
                    if(i ~= length(T))
                        if(length(T(i).imageTime) ~= T(i+1).imageData-T(i).imageData)
                            fprintf('Wrong number of images associated with Trial Cell %d Number %d\n',j,i);
                        end
                    end
                    
                    imix = ImageHelper.imageIndices(T(i));
                    if(~isempty(imix))
                        largestImageIndex = max(max(ImageHelper.imageIndices(T(i))),largestImageIndex);
                    end
                    
                    % Following are not defined for Trial 1
                    if(i > 1)
                        im1 = ImageHelper.imageIndices(T(i-1));
                        im2 = ImageHelper.imageIndices(T(i  ));
                        if(~isempty(im1) && ~isempty(im2) && im2(1) ~= im1(end)+1)
                            fprintf('Image index break between Trial Cell %d Number %d/%d (%d to %d)\n',j,i-1,i,im1(end),im2(1));
                        end
                    end
                end
                
                imageCount = zeros(largestImageIndex,1);
                for i=1:length(T)
                    imix = ImageHelper.imageIndices(T(i));
                    imageCount(imix) = imageCount(imix)+1;
                end
                if(~all(imageCount == 1))
                    fprintf(1,'Lack of a one-to-one relationship between images and indices in Trial Cell %d\n',j);
                    
                    bad = find(imageCount ==0);
                    fprintf(1,'\t%6d Images never referenced\n',numel(bad));
                    if(numel(bad) <= 50)
                        disp(bad)
                    elseif(numel(bad) > 50)
                        fprintf('Indices not printed.');
                    end
                    
                   
                    bad = find(imageCount > 1);
                    fprintf(1,'\t%6d Images repeated\n',numel(bad));
                    if(numel(bad) <= 50)
                        disp(bad)
                    elseif(numel(bad) > 50)
                        fprintf('Indices not printed.');
                    end
                end
            end
        end
        
        %% Load Images and Perform Common Processing
        
        function [S,IX,V] = convertFullToSparse(I)
            S = size(I);
            IX = find(~isnan(I));
            V = I(IX);
        end
        function I = convertSparseToFull(S,IX,V)
            I=nan(S,'single');
            I(IX) = V;
        end
        
        function [T,I,F,C] = imageLoad(filepath,opt)
            %{
                1 - Load data
                2 - Discard trials matching discard criterion (optional)
                3 - Discard images matching discard criterion (optional)
                    (Not implemented)
                4 - Convert to single precision
                5 - Convert to %-signal (optional)
                6 - Assign trials to a unique category
                7 - Regress out eye-data (optional)
                8 - Remove outlier image values (optional)
                9 - Align images to a set time-point (optional)
            %}
            if(nargin == 1)
                opt = ImageLoader.imageLoadOptions;
            end
            
            disp('Loading file.');
            T = load(filepath,'T');
            I = load(filepath,'IA');
            F = load(filepath,'F');
            
            if(~isempty(opt.trialDiscard))
                disp('Using function in opt.trialDiscard to remove trials from further analysis.');
                removeThese = false(size(T));
                for t=1:length(T)
                    if(opt.trialDiscard(T(t)))
                        removeThese(t) = true;
                    end
                end
                [T,I,~] = ImageHelper.removeTrial(T,I,[],removeThese);
            end
            
            disp('Converting to single precision with 0 values set to NaN')
            I = single(I);
            I(I<1) = nan;
            
            if(opt.percentChange)
                disp('Converting to %-signal');
                [I,T] = ImageHelper.convertToPercentChangeFromBaseline(I,T);
            end
            
            if(~isempty(opt.categorizeTrial))
                disp('Using function in opt.categorizeTrial to assign each trial a positive integer category');
                C = nan(size(T));
                for t=1:length(T)
                    try
                        C(t) = opt.categorizeTrial(T(t));
                    catch
                        % do nothing
                    end
                end
                
                % Look for invalid categories
                removeThese = find(isnan(C) | C<1 | C ~= floor(C)); 
                if(numel(removeThese) > 0)
                    disp(sprintf('%d trials returned an invalid category and were removed',numel(removeThese)));
                    [T,I,~] = ImageHelper.removeTrial(T,I,[],removeThese);
                    C(removeThese) = [];
                end
            else
                disp('Assigning all trials to the default category "1"');
                C = ones(size(T));
            end
            
            if(opt.regressEyeData)
                disp('Regressing eye data seperaely from each trial category IS NOT IMPLEMENTED');
            end
            
            if(~isempty(opt.removeOutliarImageValues))
                disp('Removing pixles with greater than opt.removeOutliarImageValues sigmas from the mean');
                mn = nanmean(I(:));
                st = nanstd(I(:));
                
                I(I < mn-opt.removeOutliarImageValues*st) = nan;
                I(I > mn+opt.removeOutliarImageValues*st) = nan;
            end
            
            
        end
        
        function opt = imageLoadOptions(varargin)
            opt.regressEyeData = true;
            
            
            for i=1:2:length(varargin)
                if(isfield(opt,varargin{i}))
                    opt.(varargin{i}) = varargin{i+1};
                else
                    error(sprintf('Invalid field "%s"',varargin{i}));
                end
            end
        end
        
        %% Remove Trials
        function [T,I,A] = removeTrial(T,I,A,i)
            if(length(i) == length(T) && islogical(i))
                fprintf(1,'Logical indexing is assumed in removeTrial, but the preferred\nmethod is explicit indexing. Logical index detection is not\ngaurenteed under all circumstances.\n');
                i = find(i);
            end
            if(isempty(T))
                return;
            end
            i = sort(i);
            if(numel(i) == 0)
                return
            end
            
            if(numel(i) == 1)
                [T,I] = ImageHelper.removeOneTrial(T,I,i);
                return
            end
            
            imA = zeros(size(T));
            imZ = zeros(size(T));
            for ii=1:numel(T)-1
                imA(ii) = T(ii).imageData;
                imZ(ii) = T(ii+1).imageData-1;
            end
            imA(end) = T(end).imageData;
            imZ(end) = size(I,3);
            
            imRemove = [];
            for ii=1:length(i)
                imRemove = [imRemove imA(i(ii)):imZ(i(ii))];
            end
            
            
            T(i) = [];
            I(:,:,imRemove) = [];
            if(~isempty(A))
                A(:,:,imRemove) = [];
            end
            imA(i) = [];
            imZ(i) = [];
            
            % Now recount images
            imC = imZ+1-imA;
            T(1).imageData = 1; % always true
            for ii=2:length(T)
                T(ii).imageData = T(ii-1).imageData+imC(ii-1);
            end
        end
        function [T,I] = removeOneTrial(T,I,i)
            
            imageIndexA = T(i).imageData;
            if(i == length(T))
                imageIndexZ = size(I,3);
            else
                imageIndexZ = T(i+1).imageData-1;
            end
            
            %fprintf('%d  %d  %d  \n',imageIndexA,imageIndexZ,size(I,3));
            
            numImagesRemoved = imageIndexZ - imageIndexA + 1;
            I(:,:,imageIndexA:imageIndexZ) = [];
            T(i) = [];
            
            for ii = i:length(T)
                T(ii).imageData = T(ii).imageData - numImagesRemoved;
            end
        end
        
        function [T,I,A] = removeTrialsShorterThan(T,I,A,time)
            ix = true(size(T));
            for i=1:length(T)
                ix(i) = ix(i) && ( ...
                    ~isempty(ImageHelper.getEndOfFixationPeriod(T(i))) ...
                    && ~isempty(T(i).stimulusOnTime) ...
                    && (ImageHelper.getEndOfFixationPeriod(T(i)) - T(i).stimulusOnTime) > time ...
                    );
            end
            [T,I,A]   = ImageHelper.removeTrial(T,I,A,find(~ix));
        end
    end
end