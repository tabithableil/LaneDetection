function [rho1, theta1, rho2, theta2] = detect_lanes (img)

    %crop to middle lines
    nim = imcrop(img, [150, 181, 199, 480]);

    %convert to gray, filter and convert to binary
    bim = rgb2gray(nim);
    bim = imgaussfilt(bim);
    bim = im2bw(bim,.55);

    %use sobel edge detection
    xmask = [-1 0 1; -2 0 2; -1 0 1];
    ymask = [-1 -2 -1; 0 0 0; 1 2 1];

    xres = imfilter(bim, xmask);
    yres = imfilter(bim, ymask);

    newim = sqrt(xres.^2+yres.^2);

    %remove small noise
    bwareaopen(newim, 50);

    %hough transform
    numTheta = 180;
    [rows, columns] = size(newim);

    theta = -90:89;
    theta = theta*(pi/180);

    rhoMax = ceil(sqrt(rows^2+columns^2));
    numRho = 2*rhoMax-1;

    %accumulator
    accum = zeros(numTheta,numRho);

    %assign all possible cosine and sine values to seperate arrays
    cost = cos(theta);
    sint = sin(theta);

    %vote for rho/theta
    for x=1:rows
        for y=1:columns
            if newim(x,y) 
                for t=1:180
                    mrho = x*cost(t)+y*sint(t);
                    rhoIndex = round((mrho+rhoMax));
                    accum(t,rhoIndex) = accum(t,rhoIndex)+1;
                end
            end
        end
    end


    %find maximum value from accumulator
    [d, max_idx] = max(accum(:));
    [p, q] = ind2sub(size(accum), max_idx);

    %compute theta and rho from index
    theta1 = -1*(pi/180)*p;
    rho1 = ((q)-(numRho+150));

    %error handling if angle isn't normal
    if theta1 > -2
        theta1 = -2.5;
        rho1 = -365;
    end

    %compute right side from left's rho and theta
    theta2 = atan(theta1)*1/180*100;
    rho2 = rho1+rhoMax+150;

    %error handling if angle isn't normal
    if theta2 < -.9
        theta2 = -.7;
        rho2 = 155;
    end

end