%function [ action ] = TrackingWebcamEaxmple
%% Object Detection and Tracking Using Live Video Acquisition
% This example shows how to automatically detect and track a Object in a live
% video stream, using the KLT algorithm.   
%
%   Copyright 2014 The MathWorks, Inc.

%% Overview
% Object detection and tracking are important in many computer vision
% applications including activity recognition, automotive safety, and
% surveillance.  In this example you will develop a simple system for
% tracking a single Object in a live video stream captured by a webcam.
% MATLAB provides webcam support through a Hardware Support Package,
% which you will need to download and install in order to run this example. 
% The support package is available via the 
% <matlab:supportPackageInstaller Support Package Installer>.
%
% The Object tracking system in this example can be in one of two modes:
% detection or tracking. In the detection mode you can use a
% |vision.CascadeObjectDetector| object to detect a Object in the current
% frame. If a Object is detected, then you must detect corner points on the 
% Object, initialize a |vision.PointTracker| object, and then switch to the 
% tracking mode. 
%
% In the tracking mode, you must track the points using the point tracker.
% As you track the points, some of them will be lost because of occlusion. 
% If the number of points being tracked falls below a threshold, that means
% that the Object is no longer being tracked. You must then switch back to the
% detection mode to try to re-acquire the Object.

%% Setup
% Create objects for detecting Objects, tracking points, acquiring and
% displaying video frames.
addpath('detectors');
addpath('generate_skinmap');
clear cam;
myObj = VideoWriter('newfile.avi');
close(myObj);
% Create the Object detector object.
faceDetector = vision.CascadeObjectDetector();
% noseDetector is more precise than face
noseDetector = vision.CascadeObjectDetector('Nose', 'UseROI', true);

% Create the point tracker object.
pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
pointTrackerThread = 10;

% Create a HOGTracker object.
HOGTracker = vision.HistogramBasedTracker;
HOGScoreThread = 0.2;
% Create the webcam object.
cam = webcam();

% Capture one frame to get its size.
videoFrame = snapshot(cam);
frameSize = size(videoFrame);

% Create the video player object. 
videoPlayer = vision.VideoPlayer('Position', [100 100 [frameSize(2), frameSize(1)]+30]);

%% Detection and Tracking
% Capture and process video frames from the webcam in a loop to detect and
% track a Object. The loop will run for 400 frames or until the video player
% window is closed.

runLoop = true;
numPts = 0;
frameCount = 0;
HOGScore = 0;

%keep the trace of the motion
handTraceQueue = CQueue;
%handtrace.time = now;     %mark the current time
%handtrace.bboxPoints = bboxPoints   %mark the center position of object
%
%record frequence
recordFrequence = 4;
countFrequence = 0;
maxQueueSize = 8;

flag = 0;

addpath('bitmap_plot_v1_2');

%write in video
myObj = VideoWriter('newfile.avi');
myObj.FrameRate = 30;
open(myObj);

scale1 = 0.1;
scale2 = 0.4;
while runLoop && frameCount < 4000
    
    % Get the next frame.
    videoFrame = snapshot(cam);
    videoFrameGray = rgb2gray(videoFrame);
    frameCount = frameCount + 1;
    countFrequence = countFrequence+1;
    %if numPts < pointTrackerThread
    if HOGScore < HOGScoreThread
        % Detection mode.
        facebox = faceDetector.step(imresize(videoFrameGray,scale2));
        if ~isempty(facebox)
            nbox = step(noseDetector, imresize(videoFrame,scale2), facebox(1,:));
        else
            nbox = [];
        end
        fbox = nbox/scale2;
        bbox = [];        
        [out bin] = generate_skinmap(videoFrame);
        %figure; imshow(bin);
        if ~isempty(fbox)
           findex1 =  fbox(2)-60;
           findex2 =  fbox(2)+fbox(4)+70;
           findex3 =  fbox(1)-60;
           findex4 =  fbox(1)+fbox(3)+50;
           if findex1 < 0
               findex1 = 1;
           end
           if findex2 > hightofFrame
               findex2 = hightofFrame;
           end
           if findex3 < 0
               findex3 = 1;
           end
           if findex4 > widthofFrame
               findex2 = widthofFrame;
           end           
           bin(findex1:findex2,findex3:findex4) =0;
        end
        L = bwlabeln(bin, 8);
        S = regionprops(L, 'Area');
        P = 1000;
        bw2 = ismember(L, find([S.Area] >= P));
        L2 = bwlabeln(bw2, 8);
        STATS = regionprops(L2,'Centroid','MajorAxisLength','MinorAxisLength');
        %figure;imshow(bw2);
        if length(STATS) >0
            maxMajorIndex = 1;
            maxMajorAxisLength = STATS(maxMajorIndex).MajorAxisLength;
            for StateIndex=1:length(STATS)
                if maxMajorAxisLength < STATS(StateIndex).MajorAxisLength
                    maxMajorIndex = StateIndex;
                    maxMajorAxisLength = STATS(StateIndex).MajorAxisLength;
                end
            end
            bbox = [floor(STATS(maxMajorIndex).Centroid(1)-(STATS(maxMajorIndex).MinorAxisLength+20)/2) floor(STATS(maxMajorIndex).Centroid(2)-(STATS(maxMajorIndex).MinorAxisLength+20)/2) floor(STATS(maxMajorIndex).MinorAxisLength+20) floor(STATS(maxMajorIndex).MinorAxisLength+20)];
        end
        
        %hand detection rely on the face detection
        if isempty(fbox)
            bbox = [];
        end
        
        if ~isempty(bbox)
%             % Find corner points inside the detected region.
%             points = detectMinEigenFeatures(videoFrameGray, 'ROI', bbox(1, :));
%             
%             % Re-initialize the point tracker.
%             xyPoints = points.Location;
%             numPts = size(xyPoints,1);
%             release(pointTracker);
%             initialize(pointTracker, xyPoints, videoFrameGray);
%             
%             % Save a copy of the points.
%             oldPoints = xyPoints;
%             
%             % Convert the rectangle represented as [x, y, w, h] into an
%             % M-by-2 matrix of [x,y] coordinates of the four corners. This
%             % is needed to be able to transform the bounding box to display
%             % the orientation of the Object.
             bboxPoints = bbox2points(bbox(1, :));  
             
             % Convert the box corners into the [x1 y1 x2 y2 x3 y3 x4 y4] 
             % format required by insertShape.
             bboxPolygon = reshape(bboxPoints', 1, []);
%             
%             % Display a bounding box around the detected Object.
             videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
%             
%             % Display detected corners.
%             videoFrame = insertMarker(videoFrame, xyPoints, '+', 'Color', 'white');
            release(HOGTracker);
            [hueChannel,~,~] = rgb2hsv(videoFrame);
            initializeObject(HOGTracker, imresize(hueChannel,scale1), ceil(scale1*abs(bbox(1,:))));
            
            % record the trace
            if countFrequence == recordFrequence
                 countFrequence = 0;
                 handtrace.time = now;
                 handtrace.bboxPoints = bboxPoints;
                 handTraceQueue.push(handtrace);
            end 
       
        end
        
    else
        % Tracking mode.
%         [xyPoints, isFound] = step(pointTracker, videoFrameGray);
%         visiblePoints = xyPoints(isFound, :);
%         oldInliers = oldPoints(isFound, :);
%                 
%         numPts = size(visiblePoints, 1);       
%         
%         if numPts >= 10
%             % Estimate the geometric transformation between the old points
%             % and the new points.
%             [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
%                 oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);            
%             
%             % Apply the transformation to the bounding box.
%             bboxPoints = transformPointsForward(xform, bboxPoints);
%             
%             % Convert the box corners into the [x1 y1 x2 y2 x3 y3 x4 y4] 
%             % format required by insertShape.
             
    [hueChannel,~,~] = rgb2hsv(videoFrame);
    % Track using the Hue channel data
    [hbox ,~, HOGScore]= step(HOGTracker, imresize(hueChannel,scale1));
    bbox = floor(hbox / scale1);
    bboxPoints = bbox2points(bbox(1, :)); 

             bboxPolygon = reshape(bboxPoints', 1, []);            
             
             % Display a bounding box around the Object being tracked.
             videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
%             
%             % Display tracked points.
%             videoFrame = insertMarker(videoFrame, visiblePoints, '+', 'Color', 'white');
%             
%             % Reset the points.
%             oldPoints = visiblePoints;
%             setPoints(pointTracker, oldPoints);

%             % record the trace
            if countFrequence == recordFrequence
                 countFrequence = 0;
                 handtrace.time = now;
                 handtrace.bboxPoints = bboxPoints;
                 handTraceQueue.push(handtrace);
            end   
%         end
  
    end
    
    action = actionRecognition(handTraceQueue);
    %if the traceQueuesize is more than maxqueuesize,pop it out
    if handTraceQueue.size() == maxQueueSize
        handTraceQueue.pop();
    end
    
    if action~=0
       handTraceQueue.empty();
       flag = action;
    end
      
    % Display the annotated video frame using the video player object.

    for fk=1:3
         finalVideoFrame(:,:,fk)=fliplr(videoFrame(:,:,fk));
    end
    
    if flag == 1
         lines={'Turn on','Hand Action 1.1'};
         finalVideoFrame=256*bitmaptext(lines,double(finalVideoFrame)/256,[1 1],struct('Color',[1 1 1 1]));
    end
    
    if flag == 2
         lines={'Turn off','Hand Action 1.1'};
         finalVideoFrame=256*bitmaptext(lines,double(finalVideoFrame)/256,[1 1],struct('Color',[1 1 1 1]));
    end
    
    if flag == 3
         lines={'left','Hand Action 1.1'};
         finalVideoFrame=256*bitmaptext(lines,double(finalVideoFrame)/256,[1 1],struct('Color',[1 1 1 1]));
    end
    
    if flag == 4
         lines={'right','Hand Action 1.1'};
         finalVideoFrame=256*bitmaptext(lines,double(finalVideoFrame)/256,[1 1],struct('Color',[1 1 1 1]));
    end
    
    step(videoPlayer, uint8(finalVideoFrame));
    writeVideo(myObj,uint8(finalVideoFrame));
    % Check whether the video player window has been closed.
    runLoop = isOpen(videoPlayer);
end

% Clean up.
clear cam;
close(myObj);
release(videoPlayer);
release(pointTracker);
release(faceDetector);

%% References
% Viola, Paul A. and Jones, Michael J. "Rapid Object Detection using a
% Boosted Cascade of Simple Features", IEEE CVPR, 2001.
%
% Bruce D. Lucas and Takeo Kanade. An Iterative Image Registration 
% Technique with an Application to Stereo Vision. 
% International Joint Conference on Artificial Intelligence, 1981.
%
% Carlo Tomasi and Takeo Kanade. Detection and Tracking of Point Features. 
% Carnegie Mellon University Technical Report CMU-CS-91-132, 1991.
%
% Jianbo Shi and Carlo Tomasi. Good Features to Track. 
% IEEE Conference on Computer Vision and Pattern Recognition, 1994.
%
% Zdenek Kalal, Krystian Mikolajczyk and Jiri Matas. Forward-Backward
% Error: Automatic Detection of Tracking Failures.
% International Conference on Pattern Recognition, 2010

displayEndOfDemoMessage(mfilename)

%end

