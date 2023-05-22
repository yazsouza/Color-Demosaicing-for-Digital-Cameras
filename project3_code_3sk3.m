% Yasmin Souza
%400325538
close all; clc;

test_image=imread('test4.png'); %input test image
%ground=imread('test3_ground.jpeg'); %ground truth image for rmse calculation
[height,length,~]=size(test_image);%test image size

%4 training images to calculate coefficient matrix
training_image1=imread('pic1.jpg');
training_image2=imread('pic2.jpg');
training_image3=imread('pic3.jpg');
training_image4=imread('pic4.jpg');

%concatenate training images
training_image = cat(2, training_image1, training_image2,training_image3,training_image4);
[height_2,length_2,~]=size(training_image); %training image size

%bayer images for test and training
bayer_image=double(bayer(test_image));
bayer_training=double(bayer(training_image));

%initiliaze values for for loop 
mosaic=zeros((height_2*length_2/4),49,4);%holds 4 matrices for 4 cases (RGGB, GRBG, GBRG, BGGR)
g_mat=zeros((height_2*length_2/4),8);%for matrix of ground truth valyes

training_s=enlarge_boarder(sum(bayer_training,3)); %make single layer and expand

%initialize incrementers for for loop
i=0;
j=0;
k=0;
l=0;

%finding mosaic patches for 4 cases (RGGB, GRBG, GBRG, BGGR)
for y =4:height_2+3 %(adding 3 to initial and final to account for expansion)
    for x =4:length_2+3
        submatrix = training_s(y-3:y+3,x-3:x+3); %7x7 submatrix
        submatrix=reshape(submatrix.',1,[]); %row
        if (mod(y,2)==0 && mod(x,2)==0)
            %RGGB
            i=i+1;
            mosaic( i, :, 1 ) = submatrix; %save to first matrix
            g_mat(i,1)=training_image(y-3,x-3,2);%ground truth matrix
            g_mat(i,2)=training_image(y-3,x-3,3);
        elseif mod(y,2)==1 && mod(x,2)==0
            %GBRG
            k=k+1;
            mosaic( k, :, 3 ) = submatrix; %save to second matrix
            g_mat(k,5)=training_image(y-3,x-3,1);
            g_mat(k,6)=training_image(y-3,x-3,3);
        elseif mod(y,2)==0 && mod(x,2)==1
            %GRBG
            j=j+1;
            mosaic( j, :, 2 ) = submatrix; %save to third matrix
            g_mat(j,3)=training_image(y-3,x-3,1);
            g_mat(j,4)=training_image(y-3,x-3,3);

        else
            %BGGR
            l=l+1;
            mosaic( l, :, 4 ) = submatrix; %save to 4th matrix
            g_mat(l,7)=training_image(y-3,x-3,1);
            g_mat(l,8)=training_image(y-3,x-3,2);
        end
    end
end


coeffs=zeros(49,8);%initialize coefficients
for i=1:8
    x_i=round(i/2);
    val=mosaic(:,:,x_i);
    x_t=transpose(val);
    coeffs(:,i)=inv(x_t*val)*x_t*g_mat(:,i); %find optimal coefficients based on training data
end

%for demosaicing 
demosaiced_image=double(bayer_image); %create double type copy of bayer image for final demosaic
bayer_single=enlarge_boarder(sum(double(bayer_image),3)); %make single layer and expand

%demosaic by X*A formula from project guide
for y =1+3:height+3
    for x =1+3:length+3
        submatrix = bayer_single(y-3:y+3,x-3:x+3); 
        submatrix=reshape(submatrix.',1,[]);
        if (mod(y,2)==0 && mod(x,2)==0)
            %RGGB
            demosaiced_image(y+-3,x+-3,2)=submatrix*coeffs(:,1); %subtracting 3 to deal with expansion
            demosaiced_image(y+-3,x+-3,3)=submatrix*coeffs(:,2);
        elseif mod(y,2)==1 && mod(x,2)==0
            %GBRG
            demosaiced_image(y+-3,x+-3,1)=submatrix*coeffs(:,5);
            demosaiced_image(y+-3,x+-3,3)=submatrix*coeffs(:,6);
        elseif mod(y,2)==0 && mod(x,2)==1
            %GRBG
            demosaiced_image(y+-3,x+-3,1)=submatrix*coeffs(:,3);
            demosaiced_image(y+-3,x+-3,3)=submatrix*coeffs(:,4);
        else
            %BGGR
            demosaiced_image(y+-3,x+-3,1)=submatrix*coeffs(:,7);
            demosaiced_image(y+-3,x+-3,2)=submatrix*coeffs(:,8);
        end
    end
end


%demosaiced image from linear regression
figure(1)
imshow(uint8(demosaiced_image));
title('Demosaiced Image Linear Regression');
%rmse = sqrt(immse(uint8(demosaiced_image),ground))%compute rmse, comment out if no ground truth

%Using matlabs demosaic function
matlab_output = demosaic(uint8(test_image), "rggb");
figure(3)
imshow(matlab_output)
title('Demosaiced Image Using Matlab')
%rmse_matlab = sqrt(immse((matlab_output),ground)) %compute MATLAB rmse


function bayer_out=bayer(image_in) %to convert to a bayer image

[height,length,~]=size(image_in); %initialize to size of image
bayer=image_in*0; %preallocate memory
for x =1:height
    for y =1:length
        if mod(x,2)==1 && mod(y,2)==1 %R
            bayer(x,y,1)=255;
        elseif mod(x,2)==0 && mod(y,2)==0 %B
            bayer(x,y,3)=255;
        else %G
            bayer(x,y,2)=255;
        end
    end
end
bayer_out= image_in.*(bayer/255);
end

function image_exp=enlarge_boarder(image)% expand pic by 3 pixels by mirroring surroundings
[height,length,~]=size(image);
image_exp=zeros(height+6,length+6);%allocate expanded space
image_exp(4:height+3,4:length+3)=image;

for y =1:3%expand top and bottom
    image_exp(y:y,1:length+6)=image_exp(8-y:8-y,1:length+6);
    image_exp(y+height+3:y+height+3,1:length+6)=image_exp(3+height-y:3+height-y,1:length+6);
end
for x =1:3%expand left and right
    image_exp(4:height+3,x:x)=image(1:height,5-x);
    image_exp(4:height+3,x+length+3:x+length+3)=image(1:height,length-x);
end

end