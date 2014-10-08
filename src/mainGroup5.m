function [ ] = mainGroup5( img )
% Reading the input image
img = imread(img);
% Converting the image to double
imdoub = im2double(img);
% Converting the double image to gray scale
imd = rgb2gray(imdoub);
figure, imshow(imd), title('Gray Level Image');

% Region Growing algorithm to identify the frames in comic page

% Taking a 5 pixel frame around the permiter of image and compute the
% average value
i = imd(:,1:5);
j = imd(:,size(imd,2)-4:size(imd,2));
k = imd(size(imd,1)-4:size(imd,1),:);
l = imd(1:5,:);
avg = (sum(sum(i)) + sum(sum(j)) + sum(sum(k)) + sum(sum(l))) / (size(imd,1)*10 + size(imd,2)*10);
avg = avg - 0.1 ;

% Selecting the four seed points at the four corners of image.
seed=false(size(imd));
seed(1,1) = true;
seed(1,size(seed,2)) = true;
seed(size(seed,1),1) = true;
seed(size(seed,1),size(seed,2)) = true;
seedf = find(seed);
seedvalues = imd(seedf);

outimg = false(size(imd));

% Region growing by comparing how close the pixel is to the avg value of 5
% pixel frame.
for itr=1:length(seedvalues)
    value =seedvalues(itr);
    final = abs(imd - value) <=  abs(avg - value);
    outimg = outimg | final;
end
figure, imshow(outimg, []), title('Binary Image after region growing');

[binimg ] = bwlabel(imreconstruct(seed, outimg),8);
% figure,imshow(binimg), title('bin img');

binimgc = imcomplement(binimg);
% figure, imshow(binimgc, []), title('complement');


% Considering the width of frame to be atleast 1/6 of image width and
% height of frame to be atleast 1/8 of image height. Multiplying both to
% get the minimum frame pixel area.
width = size(binimgc,2);
x= width/6;
x = double(floor(x));

height = size(binimgc,1);
y= height/8;
y = double(floor(y));

ar = x*y;

% All the frames having fewer than minimum area pixels will be removed
binimgc = bwareaopen(binimgc, ar);
% Removing any black holes inside the frames
binimgc = imfill(binimgc,4,'holes');
figure, imshow(binimgc) , title('After filling holes and bwareaopen');

binimgc = imcomplement(binimgc);
% Deciding the value of N i.e number of times to dilate the image followed
% by same number of erosions so as to break the link between two frames. N
% is decided based on the image width.
% Structuring element used is 3*3 containg all ones
width = size(binimgc,2);
x= width/6;
x = x*0.1;
N =int16(floor(x));

str = ones(3);

% Dilating N times
for p=1:N
    binimgc = imdilate(binimgc,str);
end
% figure, imshow(binimgc, []), title('Dilation');
% Eroding N times
for p=1:N
    binimgc = imerode(binimgc,str);
end
% figure, imshow(binimgc, []), title('Erosion');
binimgc = imcomplement(binimgc);
% figure, imshow(binimgc, []), title('After morph ');

% Using the connected component so as to separate out the frames connected
% in white.
CC = bwconncomp(binimgc);
disp(CC);
list = CC.PixelIdxList;
% Bounding box property or regionprops will be used so as to get the bounding box around the
% connected components.
s = regionprops(CC, 'BoundingBox');
disp(s);

% Loop on all the connected components(frames) and crop the components from the
% original double image with the help of bounding box
for i =1:length(s)
    fig = imcrop(imdoub, s(i).BoundingBox);
    
    %     Convert the frame extracted from rgb to hsv and continue with the
    %     value component
    hsv=rgb2hsv(fig);
    % figure,imshow(hsv(:,:,1)),title('Hue');
    % figure,imshow(hsv(:,:,2)),title('sat');
    % figure,imshow(hsv(:,:,3)),title('value');
    figv = hsv(:,:,3);
    
    %     Calculating the frame area
    g=size(fig,1);
    h=size(fig,2);
    framearea=g*h;
    
    
    %     Thresholding the Value component so as to obtain the binary image on
    %     which connected component will be applied so as to exract the
    %     balloons inside the frames.
    figt = figv > 0.9;
    % figure,imshow(figt),title('Thresholded image');
    %     Applying the connected component on the thresholded image
    cc1=bwconncomp(figt);
    disp(cc1);
    
    %     Use the bounding box , Area and Filled image property of regionprops
    s1=regionprops(cc1,'BoundingBox','Area', 'FilledImage');
    
    disp(s1);
    list = cc1.PixelIdxList;
    %     Iterate on all the connected components
    for j =1:length(s1)
        %         Get the of connected component from the Area property of bounding
        %         box and calculate the ratio of its area to frame area
        a = s1(j).Area;
        arearatio = a/framearea;
        
        %         Only those components whose ratio is greater than 0.01 are kept
        if (arearatio > 0.01)
            
            %           Get the connected component image with inner holes filled
            img = s1(j).FilledImage;
            
            %           Fill the holes if any left within the image
            img = imfill(img,'holes');
            
            ha = img(:,size(img,2));
            hb = img(:,1);
            hc = img(size(img,1),:);
            hd = img(1,:);
            
            
            %             Dilating the image and then using imfill so as to remove the
            %             black holes that are connected to the boundary of connected
            %             component
            img = imdilate(img,str);
            img = imfill(img,'holes');
            img = imdilate(img,str);
            img = imfill(img,'holes');
            img = imdilate(img,str);
            img = imfill(img,'holes');
            img = imdilate(img,str);
            img = imfill(img,'holes');
            
            img(:,size(img,2)) = ha;
            img(:,1) = hb;
            img(size(img,1),:) = hc;
            img(1,:) = hd;
            
            %             Cropping the connected component from the value image using
            %             bounding box
            fig11=imcrop(figv,s1(j).BoundingBox);
            %     figure, imshow(fig11),title('Candidate Object ');
            
            %             complementing the image followed by threshold so as to
            %             convert to binary image
            figcc = imcomplement(fig11);
            %              figure, imshow(figcc),title('Candidate Object complement');
            
            
            thresh=figcc>graythresh(figcc);
            %             figure, imshow(thresh),title('Candidate Object');
            
            % The problem with bounding box is that it may also include other parts of
            % image which are not part of bounding box as it is rectangle and connected
            % component may not always be rectangle. The image obtained from the
            % FilledImage contains the actual connected component. It is multiplied
            % with the above thresholded image so as to get rid of part which is not
            % the part of connected component.
            [n m] = size(img);
            for in=1:n
                for jm=1:m
                    thresh(in,jm) = img(in,jm).*thresh(in,jm);
                end
            end
            
            newimg = thresh;
            %           figure, imshow(thresh),title('New Candidate Object Threshold');
            
            %           There may be a lot of connected components identified which
            %           may not be the balloon. These components needs to be removed.
            %           The idea is to dilate each component so as to join the text
            %           together and the ones with larger number of white pixels
            %           compared to black are chosen as candidate balloons
            
            % To dilate the image N times, N is decided based on the component width
            % and height.
            width = size(thresh,2);
            height = size(thresh,1);
            num = max(width,height);
            x= num/6;
            
            x = x*0.15;
            N =int16(round(x));
            
            
            for p=1:N
                thresh = imdilate(thresh,str);
            end
            
            %             Getting the number of white pixels and black pixels and
            %             taking the ratio of black to white pixels.
            nWhite= sum(thresh(:));
            nBlack = numel(thresh) - nWhite;
            ratio = nBlack/nWhite;
            
            %                 figure, imshow(thresh),title('Dilation on Candidate');
            [thr thc] = size(thresh);
            %            Only those candidate balloons whose ratio of black to white
            %            pixels is less than 1 and its height greater then 15 and width
            %            greater than 50 are considered further for text extraction and
            %            the rest are discarded.
            if (ratio < 1 && thr > 15 && thc >50)
                %           figure, imshow(newimg),title('Filtered candidate objects');
                %           figure, imshow(thresh),title('imdilateText');
                
                %           Connected component is applied on the filtered candidate object
                %           so as to extract the text as connected component.
                cc2=bwconncomp(newimg);
                disp(cc2);
                s2=regionprops(cc2,'BoundingBox','Area');
                
                %           Iterate on each of the connected component
                for l=1:length(s2)
                    %       Crop the connected component from the image using
                    %       boundingbox
                    figtext=imcrop(newimg,s2(l).BoundingBox);
                    %                     figure,imshow(figtext);
                    [tr tc] = size(figtext);
                    %        Only consider those components as text whose width is
                    %        greater than 11 pixels and height is between 11 and 40
                    %        pixels
                    if (tr >11 && tc > 11 && tr < 40)
                        %     Components which satisfy the size criteria are
                        %     further filtered based on its black to white
                        %     ratio of pixels. Only those whose ratio is less
                        %     than 2.9 are considered as text component
                        tWhite= sum(figtext(:));
                        tBlack = numel(figtext) - tWhite;
                        tra = tBlack/tWhite;
                        
                        if( tra < 2.9)
                            figure,imshow(figtext),title('Text extracted ');
                        end
                    end
                end
            end
        end
    end
end
end

