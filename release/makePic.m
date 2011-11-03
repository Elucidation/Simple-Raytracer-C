clear all;
clc;
%clf;
disp('Reading from files...');
r = load('valuesr.txt');
size(r)
g = load('valuesg.txt');
size(g)
b = load('valuesb.txt');
size(b)
disp('Finished reading.');
disp('Building Picture array...');
finalPic = [];
finalPic(:,:,1) = r;
clear r;
finalPic(:,:,2) = g;
clear g;
finalPic(:,:,3) = b;
clear b;
%imshow(finalPic);
disp('Finished building array.');
disp('Writing to file...');
imwrite(finalPic,'finalPic.jpeg','JPEG');
disp('Done.');