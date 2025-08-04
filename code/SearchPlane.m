function [output, mfe,thetaf, phif] = SearchPlane(xyz)

% Input: a point cloud xyz
% Output: a vector output containing the parameter a, b, c, d of the plane
% in the form ax+by+cz+d=0
% 
try
    
    % Translation of the point cloud: barycenter in the origin
    
    trX=mean(xyz(:,1));
    trY=mean(xyz(:,2));
    trZ=mean(xyz(:,3));
    
    xyz(:,1)=xyz(:,1)-trX;
    xyz(:,2)=xyz(:,2)-trY;
    xyz(:,3)=xyz(:,3)-trZ;
    
    % Application of the Hough Transform to find the best fitting plane in
    % the Hesse form: coord contains the parameter rho, theta and phi of 
    % the equation cos(theta)sin(phi)x + sin(phi)sin(theta)y + cos(phi) +
    % - rho=0
    
    [coord, maxCoord]=PlanesHT0_search(xyz);
    
    maxCoord
    % In case, the number of fitting plane is higher than 1, the following
    % step fixes the one with lower approximation error (MFE)
    
    mAus=1;
    for j=1:size(coord,1)
        rho=coord(j,1);
        theta=coord(j,2);
        phi=coord(j,3);
        if abs(cos(phi))>10^(-3)
            x_piano = linspace(min(xyz(:,1)),max(xyz(:,1)),200);
            y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),200);
            piano = zeros(length(x_piano)*length(y_piano),3);
            for i=1:length(x_piano)
                for k=1:length(y_piano)
                    piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                    piano( (i-1)*length(x_piano)+k,2) = y_piano(k);
                    piano( (i-1)*length(x_piano)+k,3) = (rho-x_piano(i)*cos(theta)*sin(phi)-y_piano(k)*sin(phi)*sin(theta))/cos(phi);
                end
            end
        else
            if abs(cos(theta))>10^(-3)
                y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),200);
                z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),200);
                piano = zeros(length(y_piano)*length(z_piano),3);
                for i=1:length(y_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(y_piano)+k,1) = (rho-y_piano(i)*sin(theta)*sin(phi))/(sin(phi)*cos(theta));
                        piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
                        piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
                    end
                end
            else
                x_piano=linspace(min(xyz(:,1)),max(xyz(:,1)),200);
                z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),200);
                piano = zeros(length(x_piano)*length(z_piano),3);
                for i=1:length(x_piano)
                    for k=1:length(z_piano)
                        piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
                        piano( (i-1)*length(x_piano)+k,2) = (rho)/(sin(phi)*sin(theta));
                        piano( (i-1)*length(x_piano)+k,3) = z_piano(k);
                    end
                end
            end
        end
        
        [I,dist] = knnsearch(piano,xyz);
        m=MFE(xyz,dist);
        
        if m<mAus
            mAus=m;
            mfe=m;
            output=[cos(theta)*sin(phi) sin(phi)*sin(theta) cos(phi) cos(theta)*sin(phi)*trX+sin(phi)*sin(theta)*trY+cos(phi)*trZ-rho]; %cos(theta)*sin(phi)*trX+sin(phi)*sin(theta)*trY+cos(phi)*trZ
            RHO=rho;
            thetaf=theta;
            phif=phi;
%             figure
%             axis equal
%             hold on
%             scatter3(xyz(:, 1), xyz(:, 2),xyz(:, 3),'.k');scatter3(piano(:, 1), piano(:, 2),piano(:, 3),'.r');
        end
    end
    
%     n=output(1:3);
%     z=[0 0 1];
%     if norm(z-n)>5*10^(-2) && norm(z+n)>5*10^(-2)
%         R1 = rot_mat(z',n');
%         perc=0.8;
%     else
%         R1 = eye(3);
%         output(1:3)=z;
%         output(4)=trZ-rho;
%         perc=1.1;
%     end
    
%     if abs(n(3))>10^(-3)
%         x_piano = linspace(min(xyz(:,1)),max(xyz(:,1)),200);
%         y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),200);
%         piano = zeros(length(x_piano)*length(y_piano),3);
%         for i=1:length(x_piano)
%             for k=1:length(y_piano)
%                 piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
%                 piano( (i-1)*length(x_piano)+k,2) = y_piano(k);
%                 piano( (i-1)*length(x_piano)+k,3) = (RHO-x_piano(i)*n(1)-y_piano(k)*n(2))/n(3);
%             end
%         end
%     else
%         if abs(n(1))>10^(-3)
%             y_piano=linspace(min(xyz(:,2)),max(xyz(:,2)),200);
%             z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),200);
%             piano = zeros(length(y_piano)*length(z_piano),3);
%             for i=1:length(y_piano)
%                 for k=1:length(z_piano)
%                     piano( (i-1)*length(y_piano)+k,1) = (RHO-y_piano(i)*n(2))/n(1);
%                     piano( (i-1)*length(y_piano)+k,2) = y_piano(i);
%                     piano( (i-1)*length(y_piano)+k,3) = z_piano(k);
%                 end
%             end
%         else
%             x_piano=linspace(min(xyz(:,1)),max(xyz(:,1)),200);
%             z_piano=linspace(min(xyz(:,3)),max(xyz(:,3)),200);
%             piano = zeros(length(x_piano)*length(z_piano),3);
%             for i=1:length(x_piano)
%                 for k=1:length(z_piano)
%                     piano( (i-1)*length(x_piano)+k,1) = x_piano(i);
%                     piano( (i-1)*length(x_piano)+k,2) = RHO/n(2);
%                     piano( (i-1)*length(x_piano)+k,3) = z_piano(k);
%                 end
%             end
%         end
%     end
%     
%     [I,dist] = knnsearch(piano,xyz);
%     
%     xyz_aus=piano(I,:);
%     
% %     figure
% %     axis equal
% %     hold on
% %     scatter3(xyz(:, 1), xyz(:, 2),xyz(:, 3),'.k');scatter3(xyz_aus(:, 1), xyz_aus(:, 2),xyz_aus(:, 3),'.r');
% %     
%     xyz_aus=xyz_aus*R1;
%      
%     mBB = minBoundingBox(xyz_aus(:,1:2)');
%     
% %     figure
% %     axis equal
% %     hold on
% %     scatter(xyz_aus(:,1), xyz_aus(:,2),'.k');
% %     plot(mBB(1,[1:end 1]),mBB(2,[1:end 1]),'r');
%     
%     if (norm(mBB(:,1)-mBB(:,2)) < norm(mBB(:,1)-mBB(:,4)) )
%         m = (mBB(2,4)-mBB(2,1))/(mBB(1,4)-mBB(1,1));
%         alpha = atan(m);
%     else
%         m = (mBB(2,2)-mBB(2,1))/(mBB(1,2)-mBB(1,1));
%         alpha = atan(m);
%     end
%     R2 = [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
%     
%     xyz_aus(:,1:2)=xyz_aus(:,1:2)*R2;
%     
%     mBB_aus = minBoundingBox(xyz_aus(:,1:2)');
%     
%     line_mBB=[mBB_aus(1,1)*ones(200,1) linspace(mBB_aus(2,1),mBB_aus(2,2),200)'; linspace(mBB_aus(1,2),mBB_aus(1,3),200)' mBB_aus(2,2)*ones(200,1); mBB_aus(1,3)*ones(200,1) linspace(mBB_aus(2,3),mBB_aus(2,4),200)';linspace(mBB_aus(1,1),mBB_aus(1,4),200)' mBB_aus(2,1)*ones(200,1)];
%     mBB_fin=[];
%     
%     % Linea verticale sx mBB
%     [I,dist] = knnsearch(line_mBB(1:200,:),xyz_aus(:,1:2));
%     perc=0.1;
%     aus=[];
%     while size(aus,1)<1
%         aus=xyz_aus(dist<perc*mean(dist),:);
%         perc=perc+0.1;
%     end
%     
%     point1=aus(aus(:,2)==max(aus(:,2)),1:2);
%     point2=aus(aus(:,2)==min(aus(:,2)),1:2);
%     
%     if norm(mBB_aus(1:2,1)-point2')<0.1*abs(mBB_aus(2,1)-mBB_aus(2,2))
%        aus1=mBB_aus(1:2,1)';
%     else
%        aus1=point2;
%     end
%     if norm(mBB_aus(1:2,2)-point1')<0.1*abs(mBB_aus(2,1)-mBB_aus(2,2))
%        aus2=mBB_aus(1:2,2)';
%     else
%        aus2=point1;
%     end
%     if norm(aus1-aus2)<0.1*abs(mBB_aus(2,1)-mBB_aus(2,2))
%        aus1=(aus1+aus2)/2;
%        mBB_fin=[mBB_fin; aus1];
%     else
%        mBB_fin=[mBB_fin; aus1; aus2];
%     end
%     
%      % Linea orizzaontale sup mBB
%      [I,dist] = knnsearch(line_mBB(201:400,:),xyz_aus(:,1:2));
%      perc=0.1;
%      aus=[];
%      while size(aus,1)<1
%          aus=xyz_aus(dist<perc*mean(dist),:);
%          perc=perc+0.1;
%      end
%     
%     point1=aus(aus(:,1)==max(aus(:,1)),1:2);
%     point2=aus(aus(:,1)==min(aus(:,1)),1:2);
%     
%     if norm(mBB_aus(1:2,2)-point2')<0.1*abs(mBB_aus(1,2)-mBB_aus(1,3))
%        aus1=mBB_aus(1:2,2)';
%     else
%        aus1=point2;
%     end
%     if norm(mBB_aus(1:2,3)-point1')<0.1*abs(mBB_aus(1,2)-mBB_aus(1,3))
%        aus2=mBB_aus(1:2,3)';
%     else
%         aus2=point1;
%     end
%     if norm(aus1-aus2)<0.1*abs(mBB_aus(1,2)-mBB_aus(1,3))
%         if norm(mBB_fin(end,:)-aus1)>0.1*abs(mBB_aus(1,2)-mBB_aus(1,3))
%             aus1=(aus1+aus2)/2;
%             mBB_fin=[mBB_fin; aus1];
%         end
%     else
%         if norm(mBB_fin(end,:)-aus1)<0.1*abs(mBB_aus(1,2)-mBB_aus(1,3))
%             mBB_fin=[mBB_fin; aus2];
%         else
%             mBB_fin=[mBB_fin; aus1; aus2];
%         end
%     end
%     
%     % Linea verticale dx mBB
%     [I,dist] = knnsearch(line_mBB(401:600,:),xyz_aus(:,1:2));
%     perc=0.1;
%     aus=[];
%     while size(aus,1)<1
%         aus=xyz_aus(dist<perc*mean(dist),:);
%         perc=perc+0.1;
%     end   
%     
%     point1=aus(aus(:,2)==max(aus(:,2)),1:2);
%     point2=aus(aus(:,2)==min(aus(:,2)),1:2);
%     
%     if norm(mBB_aus(1:2,3)-point1')<0.1*abs(mBB_aus(2,3)-mBB_aus(2,4))
%        aus1=mBB_aus(1:2,3)';
%     else
%        aus1=point1;
%     end
%     if norm(mBB_aus(1:2,4)-point2')<0.1*abs(mBB_aus(2,3)-mBB_aus(2,4))
%        aus2=mBB_aus(1:2,4)';
%     else
%         aus2=point2;
%     end
%     if norm(aus1-aus2)<0.1*abs(mBB_aus(2,3)-mBB_aus(2,4))
%         if norm(mBB_fin(end,:)-aus1)>0.1*abs(mBB_aus(2,3)-mBB_aus(2,4))
%             aus1=(aus1+aus2)/2;
%             mBB_fin=[mBB_fin; aus1];
%         end
%     else
%         if norm(mBB_fin(end,:)-aus1)<0.1*abs(mBB_aus(2,3)-mBB_aus(2,4))
%             mBB_fin=[mBB_fin; aus2];
%         else
%             mBB_fin=[mBB_fin; aus1; aus2];
%         end
%     end
%     
%     % Linea orizzontale inf mBB
%     [I,dist] = knnsearch(line_mBB(601:800,:),xyz_aus(:,1:2));
%     perc=0.1;
%     aus=[];
%     while size(aus,1)<1
%         aus=xyz_aus(dist<perc*mean(dist),:);
%         perc=perc+0.1;
%     end
%     
%     point1=aus(aus(:,1)==max(aus(:,1)),1:2);
%     point2=aus(aus(:,1)==min(aus(:,1)),1:2);
%     
%     if norm(mBB_aus(1:2,4)-point1')<0.1*abs(mBB_aus(1,4)-mBB_aus(1,1))
%        aus1=mBB_aus(1:2,4)';
%     else
%        aus1=point1;
%     end
%     if norm(mBB_aus(1:2,1)-point2')<0.1*abs(mBB_aus(1,4)-mBB_aus(1,1))
%        aus2=mBB_aus(1:2,1)';
%     else
%         aus2=point2;
%     end
%     if norm(aus1-aus2)<0.1*abs(mBB_aus(1,4)-mBB_aus(1,1))
%         if norm(mBB_fin(end,:)-aus1)>0.1*abs(mBB_aus(1,4)-mBB_aus(1,1))
%             aus1=(aus1+aus2)/2;
%             mBB_fin=[mBB_fin; aus1];
%         end
%     else
%         if norm(mBB_fin(end,:)-aus1)<0.1*abs(mBB_aus(1,4)-mBB_aus(1,1))
%             if norm(mBB_fin(1,:)-aus2)>0.1*abs(mBB_aus(2,4)-mBB_aus(2,1))
%                 mBB_fin=[mBB_fin; aus2];
%             end
%         else
%             mBB_fin=[mBB_fin; aus1; aus2];
%         end
%     end
%     
% %     num_piano=20;
% %     x_piano = linspace(min(xyz_aus(:,1)),max(xyz_aus(:,1)),num_piano);
% %     roof = zeros(length(x_piano)^2,3);
% %     for s=1:length(x_piano)
% %         line=[x_piano(s)*ones(200,1) linspace(min(xyz_aus(:,2)),max(xyz_aus(:,2)),200)'];
% %         [I,dist] = knnsearch(line,xyz_aus(:,1:2));
% %         perc=0.1;
% %         aus=[];
% %         while size(aus,1)<1
% %             aus=xyz_aus(dist<perc*mean(dist),:);
% %             perc=perc+0.1;
% %         end
% %         roof(1+length(x_piano)*(s-1):1:length(x_piano)*s,1) = x_piano(s)*ones(length(x_piano),1);
% %         min(aus(:,2))
% %         max(aus(:,2))
% %         roof(1+length(x_piano)*(s-1):1:length(x_piano)*s,2) = linspace(min(aus(:,2)),max(aus(:,2)),num_piano);
% %     end
%     
%     mBB_fin=mBB_fin*R2';
% %     roof(:,1:2)=roof(:,1:2)*R2';
%     
%     mBB=[mBB_fin'; zeros(1,size(mBB_fin,1))];
% 
%     mBB(3,:)=mean(xyz_aus(:,3));
%     
%     mBB=R1*mBB;
% %     roof=roof*R1';
%     
%     mBB(1,:)=mBB(1,:)+trX;
%     mBB(2,:)=mBB(2,:)+trY;
%     mBB(3,:)=mBB(3,:)+trZ;
%     
%     mBB=mBB';
%     
%     piano(:,1)= piano(:,1)+trX;
%     piano(:,2)= piano(:,2)+trY;
%     piano(:,3)= piano(:,3)+trZ;
%     
%     xyz(:,1)= xyz(:,1)+trX;
%     xyz(:,2)= xyz(:,2)+trY;
%     xyz(:,3)= xyz(:,3)+trZ;
%     
%     [I,dist] = knnsearch(piano,xyz);
%     
%     roof=piano(I,:);
    
%     roof(:,1)= roof(:,1)+trX;
%     roof(:,2)= roof(:,2)+trY;
%     roof(:,3)= roof(:,3)+trZ;
    
catch
    mfe=NaN;
    output=[];
%     mBB=[];
%     roof=[];
    thetaf=[];
    phif=[];
    return
end



end