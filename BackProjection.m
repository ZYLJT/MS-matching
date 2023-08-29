function[ClearPoints_1,ClearPoints_2]=BackProjection(InteriorPoints1,InteriorPoints2,scale_value)



    N=max(InteriorPoints1(:,3));
    M=max(InteriorPoints2(:,3));
    Key_nums1=size(InteriorPoints1,1);
    Key_nums2=size(InteriorPoints2,1);
    Image_points1=zeros(Key_nums1,2);
    Image_points2=zeros(Key_nums2,2);
    for i=1:1:Key_nums1
        x=InteriorPoints1(i,1);
        y=InteriorPoints1(i,2);
        L_layer=InteriorPoints1(i,3);
        if(L_layer==1)
            x=x;
            y=y;
        else
            x=x*scale_value^(L_layer-1);
            y=y*scale_value^(L_layer-1);
        end
        Image_points1(i,:)=[x,y];
    end
    ClearPoints_1=round(Image_points1);
    for j=1:1:Key_nums2
        xx=InteriorPoints2(j,1);
        yy=InteriorPoints2(j,2);
        R_layer=InteriorPoints2(j,3);
        if(R_layer==1)
            xx=xx;
            yy=yy;
        else
            xx=xx*scale_value^(R_layer-1);
            yy=yy*scale_value^(R_layer-1);
        end
        Image_points2(j,:)=[xx,yy];
    end
    ClearPoints_2=round(Image_points2);
end