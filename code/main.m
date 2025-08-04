clc; clear all; close all

%% Facades recognition algorithm

dirname=strcat('../LAS/');
Files_cell=dir(dirname);
Files_cell(1:2)=[];
mfeFINALE=[];

for i=1:size(Files_cell,1) % for each chuck
    
    FoldNames=Files_cell(i).name;
    FoldNames=strcat(dirname, FoldNames);
    Files=dir(FoldNames);
    
    if size(Files,1)>2
        Files(1:2)=[];
        FileNames=Files(1).name;
        
        % open the .las file containing the point cloud corresponding to
        % the i-th chunk
        aus=double(string(Files_cell(i).name));
        INDEX_ID=aus;
        FileNames=strcat(FoldNames, '/', FileNames);
        lasReader = lasFileReader(FileNames);
        [ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","Classification");
        xyz=ptCloud.Location;
        labels = label2rgb(pointAttributes.Classification);
        colorData = reshape(labels,[],3);

        
        if size(xyz,1)>50           
            
            % compute the normals to each point
            dirName=strcat('build_segm/building',  int2str(aus));
            [status, msg, msgID] = mkdir(dirName);
            nameFile=strcat(dirName, '/building',  int2str(aus), '.las');
            
            
            fileID = fopen(nameFile,'w');
            fprintf(fileID,'%14.12f %14.12f %14.12f \n', xyz(1:(size(xyz,1)-1),:)');
            fprintf(fileID,'%14.12f %14.12f %14.12f \n', xyz(size(xyz,1),:)');
            fclose(fileID);
            
            outFile=strcat(dirName, '/building_norm',  int2str(aus), '.xyz');
            command = strcat('cgal-scaline-normals.exe -i',{' '}, FileNames,' -o ',{' '}, outFile);
            command=command{1,1};
            status = system( sprintf(command) );
            
            outFile=strcat(dirName, '/building_norm',  int2str(aus), '_angle_and_flag.xyz');
            [x, y, z]=xyzread(outFile);
            nx=x(2:2:end);
            x=x(1:2:end);
            ny=y(2:2:end);
            y=y(1:2:end);
            nz=z(2:2:end);
            z=z(1:2:end);
            
            xyz=[x y z];
            normals=[nx ny nz];
            
            trX=mean(xyz(:,1));
            trY=mean(xyz(:,2));
            trZ=mean(xyz(:,3));
            xyz(:,1)=xyz(:,1)-trX;
            xyz(:,2)=xyz(:,2)-trY;
            xyz(:,3)=xyz(:,3)-trZ;
            

            xyz_orig=xyz;
           
            % Converting .las into .txt -- > x y z nx ny nz
            xyz_norm=[xyz normals];
            nameFile=strcat(dirName, '/building',  int2str(aus), '.txt');
            fileID = fopen(nameFile,'w');
            fprintf(fileID,'%14.12f %14.12f %14.12f %14.12f %14.12f %14.12f \n', xyz_norm(1:(size(xyz_norm,1)-1),:)');
            fprintf(fileID,'%14.12f %14.12f %14.12f %14.12f %14.12f %14.12f \n', xyz_norm(size(xyz_norm,1),:)');
            fclose(fileID);
            
            dirName=strcat('../output/building',  int2str(aus));
             [status, msg, msgID] = mkdir(dirName);
            outFile=strcat(dirName, '/building_output',  int2str(aus), '.ply');
            num=size(xyz,1)/5000;
            
            if num<10
                num=10;
            end 
            
            % segmentation of the point cloud into planes
            command = strcat('fit_region_growing-fit.exe -i',{' '}, nameFile,' -o ',{' '}, outFile, ' -a 10 -s',{' '}, int2str(num));
            command=command{1,1};
            status = system( sprintf(command) );
                      
            obj=pcread(outFile);
            
            xyz=obj.Location;
            color=obj.Color;
            
            mBB_max_X=max(xyz(:,1));
            mBB_min_X=min(xyz(:,1));
            mBB_max_Y=max(xyz(:,2));
            mBB_min_Y=min(xyz(:,2));
            
            control=double(color);
            control=unique(control,'rows');
            
            segment={};
            aus=[0 0 0];
            for s=1:size(control,1)
                if norm(aus-control(s,:))>0
                    segment=[segment; xyz(find(sum(color==control(s,:),2)>2),:)];
                    aus=double(control(s,:));
                end
            end

            % HT-based recognition algorithm: to each segment, an array is
            % associated: [a b c d mfe theta phi], where:
            % - a, b, c and d are the coefficient of the implicit plane
            % equation;
            % - mfe is the mean fitting error;
            % - theta and phi are the inclination angles of the plane
            segmentation={};
            Final=zeros(size(segment,1),4);
            angles=zeros(size(segment,1),2);
            mfe_final=zeros(size(segment,1),1);
            lab_vert=[];
            angles_orig=zeros(size(segment,1),2);
            facade=[];
            for s=1:size(segment,1)
                segm_xyz=segment{s,1};
                if size(segm_xyz,1)>100
                    if size(segm_xyz,1)<500
                        segm_xyz=segm_xyz(1:5:end,:);
                    else
                        segm_xyz=segm_xyz(1:10:end,:);
                    end
                end
                [output, mfe, theta, phi]=SearchPlane(segm_xyz);
                Final(s,:)=output;
                mfe_final(s)=mfe;
                angles(s,1)=theta*180/pi;
                angles(s,2)=phi*180/pi;
                angles_orig(s,1)=theta*180/pi;
                angles_orig(s,2)=phi*180/pi;

                % with opportune transformation, the angle phi correpsonds
                % to the inclination of the plane with respect the
                % horizontal plane
                if angles(s,2)>90
                    angles(s,2)=180-angles(s,2);
                    if  angles(s,1)<180
                        angles(s,1)=angles(s,1)+180;
                    else
                        angles(s,1)=angles(s,1)-180;
                    end
                end
                % in case the angle phi is > of 85 degrees, the plane is
                % vertical
                if angles(s,2)>85
                    facade=[facade;segment{s,1}];
                    segmentation=[segmentation; segment{s,1}];
                end

            end
            
            % a .txt file that cointains the label associated to each point
            % of the original point cloud is created: 1 if the point
            % belongs to a facade, 0 otherwise
            label=zeros(size(xyz_orig,1),1);
            if size(facade,1)>0
                [I,D] = knnsearch(xyz_orig,facade);
                label(I(D<10^(-3)))=1;
            end

            nameFile=strcat(dirName, '/label.txt');
            fileID = fopen(nameFile,'w');
            writematrix(label,nameFile,'Delimiter','tab')
            fclose(fileID);
            
        end
    end
    
end

