% Author: Chew Lik Siang
% Date: 2022/02/13 (last modified)
% Revision: 1
% Function: 1) To extract resonance frequency readings from screen-recorded video. 
%           2) To extract relative humidity readings (7-segments) from recorded video.
%           3) Write the data points into an Excel file.
% Procedures: 1) Click and drag to crop the area of the RF reading.
%                Be careful not to crop anything other than the interested
%                readings as that might affect the recognition. 
%             2) Click and drag to crop the area of the TEMPERATURE
%                reading(only the numerical part).
%             3) Select the segments by clicking once on each segment of
%                the digit. eg: (a1,b1...,g1,a2,b2..,g2,a3,b3..etc)
%                Recognition process will start automatically after all 3
%                segments are entered.
%             4) 3 plots will be given at the end of the scripts. 
%             5) Extracted data points are tabulated in the respective
%                Excel file.

%% User Initialisation
clc; clear all; close all;
obj_na = VideoReader('Sample 6 - NA.mp4'); % network analyser video
obj_humid = VideoReader('Sample 6 - RH.mp4'); % humidity video
weird_zero = 0; % set to 1 if the 0 printing is very similar to 9 for temperature
humid_movement_tracking = 0; % set it to 1 if the humidity reading is moving
seven_segment_threshold = 0.1; % should be < 1
video_extraction_interval = 20; 
% interval between each extraction point (bigger --> shorter elapsed time)
% e.g. video_extraction_interval = 1 --> 1 second interval

%% frequency 
iter = 1;
this_frame = read(obj_na,iter);
figure; imshow(this_frame);
roi = round(getPosition(imrect));
Icropped = imcrop(this_frame,roi);
figure; imshow(Icropped);

%% temperature
iter_temperature = 1;
this_frame_temperature = read(obj_humid,iter_temperature);
figure; imshow(this_frame_temperature);
roi_temperature = round(getPosition(imrect));
Icropped_temperature = imcrop(this_frame_temperature,roi_temperature);
figure; imshow(Icropped_temperature);

%% humidity_measurement
num_dig=3; %number of digits in the image
ntotal=num_dig*7; %total segment number
k = 1;
this_frame_humid = read(obj_humid,k);
imshow(this_frame_humid)

point = round(ginput(ntotal));%selet segments(a1,b1...,g1,a2,b2..,g2,a3,b3..etc)
x=point(:,2); 
y=point(:,1);

resize_scale = 20;
%% ROI the word "HUMIDIFY"
if humid_movement_tracking == 1
    figure; imshow(this_frame_humid);
    roi_humidify = round(getPosition(imrect));
    Icropped_humid = imcrop(this_frame_humid, roi_humidify);
    figure; imshow(Icropped_humid);

    Icropped_resize = imresize(Icropped_humid,resize_scale);
    ocrResults   = ocr(Icropped_resize);
    Iocr         = insertObjectAnnotation(Icropped_resize, 'rectangle', ...
                   ocrResults.WordBoundingBoxes, ...
                   ocrResults.WordConfidences);
    figure; imshow(Iocr);

    expected_length = length(1:round(obj_humid.FrameRate)*video_extraction_interval:obj_humid.NumFrames);
    reference_x = zeros(1, expected_length);
    reference_y = zeros(1, expected_length);
    reference_x_change = zeros(1, expected_length);
    reference_y_change = zeros(1, expected_length);
    reference_x(1) = round(ocrResults.WordBoundingBoxes(1,1)/resize_scale);
    reference_y(1) = round(ocrResults.WordBoundingBoxes(1,2)/resize_scale);
    reference_x_change(1) = 0;
    reference_y_change(1) = 0;
end

%% Enlarge for better recognition
% frequency
Icropped_resize = imresize(Icropped,5);
figure; imshow(Icropped_resize);

% temperature
Icropped_resize_temperature = imresize(Icropped_temperature,5);
figure; imshow(Icropped_resize_temperature);

%% recognition
% frequency
ocrResults = ocr(Icropped_resize,'TextLayout', 'Block', 'CharacterSet', '-0123456789.');
value = str2double(ocrResults.Text(~isspace(ocrResults.Text)));
fprintf('%d\n',value)

% temperature
ocrResults_temperature = ocr(Icropped_resize_temperature,'TextLayout', 'Block', 'CharacterSet', '-0123456789.');
value_temperature = str2double(ocrResults_temperature.Text(~isspace(ocrResults_temperature.Text)));
fprintf('%d\n',value_temperature)

%% Extract frequency readings at 1 second interval of the video
tic
iter = 1;
expected_length = length(1:round(obj_na.FrameRate)*video_extraction_interval:obj_na.NumFrames);
rf_reading = zeros(1, expected_length);
fprintf("\n")
h = waitbar(0,'Extracting Resonance Frequency...');
for img = 1:round(obj_na.FrameRate)*video_extraction_interval:obj_na.NumFrames
    this_frame = read(obj_na,img);
    Icropped = imcrop(this_frame,roi);
    Icropped_resize = imresize(Icropped,5);
    ocrResults = ocr(Icropped_resize,'TextLayout', 'Block', 'CharacterSet', '-0123456789.');
    value = str2double(ocrResults.Text(~isspace(ocrResults.Text)));
    rf_reading(iter) = value;
    iter = iter + 1;
    waitbar(img/obj_na.NumFrames,h);
end
delete(h)
toc
%% Extract temperature readings at 1 second interval of the video
tic
iter_temperature = 1;
expected_length = length(1:round(obj_humid.FrameRate)*video_extraction_interval:obj_humid.NumFrames);
temperature_reading = zeros(1, expected_length);
fprintf("\n")
h = waitbar(0,'Extracting temperature...');
SE = strel('square',4);
SE2 = strel('disk',4);
for img_temperature = 1:round(obj_humid.FrameRate)*video_extraction_interval:obj_humid.NumFrames
    this_frame_temperature = read(obj_humid,img_temperature);
    Icropped_temperature = imcrop(this_frame_temperature,roi_temperature);
    Icropped_resize_temperature = imresize(Icropped_temperature, 5);
    
    if weird_zero == 1
        Icropped_resize_temperature = bw_dilation_erosion_func(Icropped_resize_temperature);
    end

    ocrResults_temperature = ocr(Icropped_resize_temperature,'TextLayout', 'Block', 'CharacterSet', '-0123456789.');
    value_temperature = str2double(ocrResults_temperature.Text(~isspace(ocrResults_temperature.Text)));

    temperature_reading(iter_temperature) = value_temperature;
    iter_temperature = iter_temperature + 1;
    waitbar(img_temperature/obj_humid.NumFrames,h);
end
delete(h)
toc

%% threshold computation - differentiate between bright & dark segments
threshold = 0;
for j = 1:num_dig
    for i = 1:7
        pixel_value_1st(i) = this_frame_humid(x(i+(j-1)*7),y(i+(j-1)*7));
    end
    max_pixel = max(pixel_value_1st);
    min_pixel = min(pixel_value_1st);
    threshold_temp = (max_pixel-min_pixel)*seven_segment_threshold + min_pixel;
    
    if (threshold_temp - min_pixel)>10 && (threshold_temp>threshold)
        %% to get highest threshold & prevent all pixels lit up situation
        threshold = threshold_temp;
    end
end


%% Extract R.humidity readings at 1 second interval of the video
tic
expected_length = length(1:round(obj_humid.FrameRate)*video_extraction_interval:obj_humid.NumFrames);
rh_reading = zeros(1, expected_length);
pixel_value = zeros(1, 7);
whole_digit = zeros(1, 3);
h = waitbar(0,'Extracting RH...');
iter = 1;
for img = 1:round(obj_humid.FrameRate)*video_extraction_interval:obj_humid.NumFrames
    this_frame_humid = read(obj_humid, img);
    
    for j = 1:num_dig
        for i = 1:7
            pixel_value(i) = this_frame_humid(x(i+(j-1)*7),y(i+(j-1)*7));
        end

        pixel_value(pixel_value <= threshold) = 0;
        pixel_value(pixel_value >= threshold) = 1;
        digit_binary = num2str(pixel_value);
        digit_binary = digit_binary(~isspace(digit_binary));
        digit_saver = bin2dec(digit_binary);

        if (digit_saver == 126)
            whole_digit(j) = '0';

        elseif (digit_saver == 48)
            whole_digit(j) = '1';

        elseif (digit_saver == 109)
            whole_digit(j) = '2';

        elseif (digit_saver == 121)
            whole_digit(j) = '3';

        elseif (digit_saver == 51)
            whole_digit(j) = '4';

        elseif (digit_saver == 91)
            whole_digit(j) = '5';

        elseif (digit_saver == 95)
            whole_digit(j) = '6';

        elseif (digit_saver == 112)
            whole_digit(j) = '7';

        elseif (digit_saver == 127)
            whole_digit(j) = '8';

        elseif (digit_saver == 123)
            whole_digit(j) = '9';
        end
    end

    new_whole_digit = [whole_digit(1:end-1) '.' whole_digit(end)];
    digit_str = str2double(strcat(new_whole_digit));
    rh_reading(iter) = digit_str;

    if humid_movement_tracking == 1
        %% compute new reference point
        Icropped_humid = imcrop(this_frame_humid, roi_humidify);
        Icropped_resize_humid = imresize(Icropped_humid,resize_scale);
        ocrResults   = ocr(Icropped_resize_humid);

        try
            reference_x(iter+1) = round(ocrResults.WordBoundingBoxes(1,1)/resize_scale);
            reference_y(iter+1) = round(ocrResults.WordBoundingBoxes(1,2)/resize_scale);
        catch
            reference_x(iter+1) = reference_x(iter);
            reference_y(iter+1) = reference_y(iter);
            warning('OCR no result');
        end
        reference_x_change(iter) = reference_x(iter+1) - reference_x(iter);
        reference_y_change(iter) = reference_y(iter+1) - reference_y(iter);
        if abs(reference_x_change(iter))>5
            reference_x_change(iter) = 0;
        end

        if abs(reference_y_change(iter))>5
            reference_y_change(iter) = 0;
        end

    %     y = y + reference_x_change(iter);
        x = x + reference_y_change(iter);
    end

    iter = iter + 1;
    waitbar(img/obj_humid.NumFrames,h);
end
delete(h)
toc

%% write to excel
rh_reading_length = length(rh_reading);
rf_reading_length = length(rf_reading);
temperature_reading_length = length(temperature_reading);
minimum_length = min([rh_reading_length, rf_reading_length, temperature_reading_length]);
time = 1:minimum_length;
time = (time/60)*video_extraction_interval;
T = table(time', rf_reading(1:minimum_length)', rh_reading(1:minimum_length)', temperature_reading(1:minimum_length)');
T.Properties.VariableNames = {'Time (minute)', 'Resonance Frequency (MHz)', 'Relative Humidity (%)', 'Temperature (°C)'};
excel_filename = sprintf('Data_%s.xlsx', datestr(now,'mm-dd-yyyy HH-MM')); % name of the output excel file
writetable(T, excel_filename);

%% plotting
figure
subplot(3,1,1);
time = 1:length(rh_reading);
time = (time/60)*video_extraction_interval;
plot(time,rh_reading)
title("Relative humidity (%)");
xlabel("Time (min)");

subplot(3,1,2);
time = 1:length(rf_reading);
time = (time/60)*video_extraction_interval;
plot(time,rf_reading)
title("Resonance frequency (MHz)");
xlabel("Time (min)");

subplot(3,1,3);
time = 1:length(temperature_reading);
time = (time/60)*video_extraction_interval;
plot(time,temperature_reading)
title("Temperature (°C)");
xlabel("Time (min)");

rh_reading_length = length(rh_reading);
rf_reading_length = length(rf_reading);
minimum_length = min(rh_reading_length, rf_reading_length);
figure
plot(rf_reading(1:minimum_length), rh_reading(1:minimum_length))
title("Relative humidity (%) vs Resonance frequency (MHz)");
ylabel("Relative humidity (%)");
xlabel("Resonance frequency (MHz)");

figure
plot(rh_reading(1:minimum_length), rf_reading(1:minimum_length))
title("Resonance frequency (MHz) vs Relative humidity (%)");
ylabel("Resonance frequency (MHz)");
xlabel("Relative humidity (%)");

%% functions to be called in body
function output = bw_dilation_erosion_func(resized_img)
    SE = strel('square',4);
    SE2 = strel('disk',4);
    Icropped_resize_temperature = im2bw(resized_img);
    Icropped_resize_temperature = imdilate(Icropped_resize_temperature, SE);
    Icropped_resize_temperature = imerode(Icropped_resize_temperature, SE2);
    output = Icropped_resize_temperature;
end