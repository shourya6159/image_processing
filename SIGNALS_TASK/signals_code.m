% LOAD CORRECT FOLDER PATH TO READ IMAGE:
folder_path=pwd
imds = imageDatastore(fullfile(folder_path, "cursed_schematic_0*.png"));


% Loop until the datastore is out of images:
while hasdata(imds)
    [img, info] = read(imds);
    img = im2double(img);
    img=im2gray(img);
    size(img);
    [~, base_name, ~] =fileparts(info.Filename); 
    file_num =extractAfter(base_name, "cursed_schematic_");
    fprintf('Processing %s\n', base_name);
    
    
    %Removing salt and pepper noise:
    %Applying a median filter for removing salt and pepper noise
    img_filtered=medfilt2(img,[3 3]);
    
    

    
    
    % SUPPRESSING NOISE IN FREQUENCY DOMAIN:
    
    %Converting image to frequency domain:
    img_filtered_freq=fftshift(fft2(double(img_filtered)));
    [m,n]=size(img_filtered_freq);
    
    % INSIGHTS:
    % The grid noise creates a + in the frequency domain which we must remove
    % The nature of the + is such that it would pass through the center
    % The only thing we should tweak is the width of the +
    middle_row=floor(m/2)+1;
    middle_column=floor(n/2)+1;
    
    %Defining a checkpoint index that's neither too far from center nor too close
    checkpoint=128+50;
    
    %Look at the array of pixel values at checkpoint row and column:
    row_checkpoint=img_filtered_freq(checkpoint,128:256);
    column_checkpoint=img_filtered_freq(128:256,checkpoint)';
    
    % We can be sure that the brightest noise pixels will be approximately at
    % the center of the image!
    
    % Find the brightest pixel value at checkpoint row:
    brightest_row_noise_mag=abs(img_filtered_freq(checkpoint,128));
    
    % Find the brightest pixel value at checkpoint column:
    brightest_column_noise_mag=abs(img_filtered_freq(128,checkpoint));
    
    % Now lets see how far from the center will brightness continue:
    
    
    % For finding width of the + , I defined a find_threshold function at bottom
    % Scroll down to see!
    
    % AUTOMATED PEAK DETECTION BASED ON INSIGHTS: 
    % Find the width of + along x axis:
    plus_width_along_h=find_threshold(row_checkpoint,brightest_row_noise_mag,0.7);
    
    % Find the width of + along y axis:
    plus_width_along_v=find_threshold(column_checkpoint,brightest_column_noise_mag,0.7);
    
    % If we got rid of the + entirely then it would darken the entire image.
    % Since + passes through center, removing it would also remove the bright
    % DC component at center which governs the average brightness of entire
    % image
    % So I kept the central part of the + as it is by defining a DC threshold:
    DC_threshold=20;
    
    %Creating the frequency mask: (1 everywhere but dark where the + lies)
    mask=ones(m,n);
    darkening_factor=0.3 % We darken the noise by a multiplier of 0.3
    mask(floor(m/2)+DC_threshold:n,middle_column-plus_width_along_h:middle_column+plus_width_along_h)=darkening_factor;
    mask(1:floor(m/2)-DC_threshold,middle_column-plus_width_along_h:middle_column+plus_width_along_h)=darkening_factor;
    mask(middle_row-plus_width_along_v:middle_row+plus_width_along_v,1:floor(n/2)-DC_threshold)=darkening_factor;
    mask(middle_row-plus_width_along_v:middle_row+plus_width_along_v,floor(n/2)+DC_threshold:n)=darkening_factor;
    
    % Apply the mask by element-wise multiplication
    filtered_freq=img_filtered_freq.*mask;
    
    % Take the inverse fft to get back into spatial domain
    cleaned_img=ifft2(ifftshift(filtered_freq));
    cleaned_img=abs(cleaned_img);
    
    
    
    
    % APPLYING MORPHOLOGICAL OPEARTIONS TO REMOVE GRID COMPLETELY:
    
    
    % Although the grid noise has been suppressed, artifacts still remain
    % I try to remove grid completely by closing -> opening -> closing 
    % (read report to understand why)
    % The size of the structuring element is set to be large enough to eat up tiny
    % black pixels and white pixels of the grid
    grid_frequency=find_grid_frequency(img_filtered)
    grid_size=find_grid_size(grid_frequency)
    pixel =grid_size             % Not hard coded :>
    se =strel('disk', pixel);  
    
    %Firstly, close the image, the grid artifacts that persist as tiny black
    %pixels get removed
    cleaned_img_old=cleaned_img;
    cleaned_img=imclose(cleaned_img,se);
    
    %Secondly, open the image, the grid artifacts that persist as tiny white
    %pixels get removed
    %Since our object is white, I set the maximum size of the se to 4 so that opening does not eat any part of the object itself
    se_smaller=strel('disk',min(4,pixel)); 
    cleaned_img=imopen(cleaned_img,se_smaller);
    
    %Finally, close the image again with a bigger se, to remove the larger black pixels formed after closing -> opening 
    se=strel('disk',pixel+1); 
    cleaned_img=imclose(cleaned_img,se);
    
    
    
    % Displaying 5 panel output
    figure("Name", "Restored " + file_num, "NumberTitle", "off");
    subplot(2,3,1);
    imshow(img,[]);
    title('Original image');
    
    subplot(2,3,2);
    imshow(img_filtered,[]);
    title('Median filtered');
    
    subplot(2,3,3);
    imshow(log(1 + abs(img_filtered_freq)),[]);
    title('FFT spectrum');
    
    subplot(2,3,4);
    imshow(mask,[]);
    title('Frequency mask');
    
    subplot(2,3,5);
    imshow(cleaned_img,[]);
    title('Restored image');
    
    
    % Bonus: Convert image to blueprint style
    figure("Name","Blueprint "+file_num,"NumberTitle","off");
    blueprint=generate_blueprint(cleaned_img);
    imshow(blueprint,[]);
    title("Blueprint");
    drawnow;
end

%FUNCTIONS USED: 

function threshold_width = find_threshold(checkpoint, brightest_noise_mag, sensitivity)
% Finds the index where magnitude drops below a threshold
% checkpoint: The array of FFT magnitudes
% brightest_noise_mag: The reference peak magnitude
% sensitivity: The multiplier

    threshold_width = 1; %Default if no threshold is met
    
    % If the brightness drops significantly that is our threshold for the mask!
    % Traverse the row_checkpoint or column_checkpoint:
    for i = 1:length(checkpoint)
        current_pixel_mag = abs(checkpoint(i));
    
        % Update peak if we find a higher one:
        if current_pixel_mag > brightest_noise_mag
            brightest_noise_mag = current_pixel_mag;
        end
        % Check against the threshold
        if current_pixel_mag < sensitivity * brightest_noise_mag
            threshold_width = i;
            break; % Exit once the threshold is crossed
        end
end
end




function blueprint = generate_blueprint(cleaned_img)
% Converts the cleaned image into blueprint format
% Extract edges from cleaned image
    edge_mask = edge(cleaned_img, 'canny',0.125);
    
    
    [rows, cols] = size(cleaned_img);
    
    % Initializing the 3 color channels (Red, Green, Blue) with the blueprint color
    % Using (R,G,B)=(0.1,0.3,0.6) as the blueprint color
    R= zeros(rows,cols)+0.1;
    G= zeros(rows,cols)+0.3;
    B= zeros(rows,cols)+0.6;
    
    %Dilated the edges so that lighter parts which are closely spaced can stay
    %connected
    se= strel('square', 7); 
    edge_mask= imdilate(edge_mask, se);
    
    %Slightly eroded the edges so that the edges are not too thick
    se = strel('square', 5); 
    edge_mask= imerode(edge_mask,se);
    
    
    %Wherever edge_mask is true (1), set the RGB values to 0.6 (Whitish)
    R(edge_mask) = 0.6;
    G(edge_mask) = 0.6;
    B(edge_mask) = 0.6;
    
    %Concatenate the channels to form the final RGB image
    blueprint = cat(3, R, G, B);
    
end


function grid_size= find_grid_size(grid_frequency)
%Takes in the frequency of grid lines in an image and outputs the 'scaled' size of
%the grid to be used for determining the radius of the structuring element
    %I found the conditions experimentally.
    if (grid_frequency>36)
        grid_size=3;  %high grid frequency small grid_size
    end
    if (grid_frequency>=25) && (grid_frequency<36) 
        grid_size=4;
    end
    if (grid_frequency>=23) && (grid_frequency<25)
        grid_size=5;
    end
    if (grid_frequency<=22) && (grid_frequency>=17) 
        grid_size=6;
    end
    if grid_frequency<17
        grid_size=7;  %low grid frequency large grid_size
    end
end


function grid_frequency = find_grid_frequency(img_filtered)
%Takes in the median filtered image and finds the frequency of grid lines
%in it at row 4.
    row_4_pixels=img_filtered(4,:);
    black_threshold=0.23;
    low_pixels = row_4_pixels < black_threshold;
    A=low_pixels;
    num_black=sum(A(2:end)==1 & A(1:end-1)==0)
    grid_frequency=num_black;
end
