
eval( [ 'load index;' ] );
eval(['root_directory_contents = dir;']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need modifiy %
foldername1 = 'BeforeTile';
foldername2 = 'AfterTile';
ch_num =4;
%%%% create 'xyz' variable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check=[]; c=1;

numberOfContents = size(root_directory_contents,1);
numberOfDirectory = 0;
orderOfDirectory = [];
for x=1:numberOfContents
    if getfield(root_directory_contents(x), 'isdir') == 1
        numberOfDirectory = numberOfDirectory + 1;
        orderOfDirectory(numberOfDirectory) = x;
    end
end
extract_Directory = 1:numberOfDirectory-2;


mkdir(foldername1);
tile_num=zeros(500,1);

for n = extract_Directory
    
    currentDirectory = getfield(root_directory_contents(orderOfDirectory(n+2)), 'name')
    if n + 2 < numberOfDirectory
    nextDirectory = getfield(root_directory_contents(orderOfDirectory(n+3)), 'name')
    end
    
    display(currentDirectory);
    eval( [ 'cd ''' currentDirectory ''''] );
%     cd('Average10');
    
    eval(['directory_contents = dir;']);
    numberOfContents = size(directory_contents,1);
    numberOfFiles = 0;
    orderOfFiles = [];

for x=1:numberOfContents
    if getfield(directory_contents(x), 'isdir') == 0
        numberOfFiles = numberOfFiles + 1;
        orderOfFiles(numberOfFiles) = x;
    end
end

numberOfFiles=numberOfFiles-3; %ThorImage

for x=1:numberOfFiles

fname = directory_contents(orderOfFiles(x)).name;
display(fname);
imageType = fname(end-2:end);  
current_image = imread(fname,imageType);
image_x=max(xyz(:,1))+size(current_image,2); %x size [pixels]
image_y=max(xyz(:,2))+size(current_image,1); %y size [pixels]

    if x <= numberOfFiles/ch_num, z=x+xyz(n,3); ch=1; tile_num(z)=tile_num(z)+1;
    elseif x <= 2*numberOfFiles/ch_num, z=x+xyz(n,3)-numberOfFiles/ch_num; ch=2;
    elseif x <= 3*numberOfFiles/ch_num, z=x+xyz(n,3)-2*numberOfFiles/ch_num; ch=3; 
    else z=x+xyz(n,3)-3*numberOfFiles/ch_num; ch=4; 
    end 
    
    check(c,:)=[n, x, ch, z ]; c=c+1;
%     if z>z_max, z_max=z; end
    
    canvas=zeros(image_y, image_x);
    canvas(xyz(n,2)+1:xyz(n,2)+size(current_image,1),xyz(n,1)+1:xyz(n,1)+size(current_image,2))=current_image;
  
    
%     cd ..; 
    cd ..; cd(foldername1);
    imagename=['BeforeTile-' num2str(ch) '-' index(z,:) '-' index(n+1,:) '.' imageType];
    image=uint16(canvas);
    imwrite(image,imagename,imageType);
    cd ..; eval( [ 'cd ''' currentDirectory ''''] ); 
%     cd('Average10');

end
% cd ..; 
cd ..;
end


cd(foldername1);
eval(['directory_contents = dir;']);
    numberOfContents = size(directory_contents,1);
    numberOfFiles = 0;
    orderOfFiles = [];

for x=1:numberOfContents
    if getfield(directory_contents(x), 'isdir') == 0
        numberOfFiles = numberOfFiles + 1;
        orderOfFiles(numberOfFiles) = x;
    end
end

mkdir(foldername2);
z=1; y=1;
for x=1:numberOfFiles
    fname = directory_contents(orderOfFiles(x)).name;
    display(fname);
    imageType = fname(end-2:end);  
    current_image = imread(fname,imageType);
    
    if x <= numberOfFiles/ch_num, y=x; % ch=1 
    elseif x <= 2*numberOfFiles/ch_num, y=x-numberOfFiles/ch_num; % ch=2 
    elseif x <= 3*numberOfFiles/ch_num, y=x-2*numberOfFiles/ch_num; % ch=3 
    else y=x-3*numberOfFiles/ch_num; % ch=4 
    end 
    
    if y==1, z=1;end

    if y==1 || y==sum(tile_num(1:z-1))+1, image=current_image;
    elseif y<=sum(tile_num(1:z)), image=max(image, current_image);
    end
    if y==sum(tile_num(1:z)) 
        cd(foldername2); 
        imagename=[fname(12:end-4) '.' imageType]; 
        image16=uint16(image);
        imwrite(image16,imagename,imageType);
        cd ..; 
        z=z+1;
    end    
    
end
cd ..;

