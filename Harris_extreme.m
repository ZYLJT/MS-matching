function key_point_array=Harris_extreme(Harris_function,gradient,angle,layer,threshold,sigma,ratio)


    if nargin<4
        layer=3;
    end
    if nargin<5
        threshold=0.4;
    end
    if nargin<6
        sigma=1.6;
    end
    if nargin<7
        ratio=2^(1/3);
    end



    BORDER_WIDTH=4;
    HIST_BIN=36;
    HAPCG_ORI_PEAK_RATIO=0.8;
    key_number=0;
    %first_points=1500;
    first_points=100;
    key_num=layer*first_points;

    key_point_array=zeros(key_num,4);
    for i=1:1:layer
        temp_current=Harris_function{i};
        gradient_current=gradient{i};
        angle_current=angle{i};


        M=size(temp_current,1);
        N=size(temp_current,2);
        for j=BORDER_WIDTH:1:M-BORDER_WIDTH
            for k=BORDER_WIDTH:1:N-BORDER_WIDTH
                temp=temp_current(j,k);
                if(temp>(threshold*1.4)&&...
                    temp>temp_current(j-1,k-1)&&temp>temp_current(j-1,k)&&temp>temp_current(j-1,k+1)&&...
                    temp>temp_current(j,k-1)&&temp>temp_current(j,k+1)&&...
                    temp>temp_current(j+1,k-1)&&temp>temp_current(j+1,k)&&temp>temp_current(j+1,k+1))
                    scale=sigma*ratio^(i-1);

                    [hist,max_value]=Hist_Oritation(k,j,scale,gradient_current,angle_current,HIST_BIN);
                    mag_thr=max_value*HAPCG_ORI_PEAK_RATIO;
                    for kk=1:1:HIST_BIN
                        if(kk==1)
                            k1=HIST_BIN;
                        else
                            k1=kk-1;
                        end
                        if(kk==HIST_BIN)
                            k2=1;
                        else
                            k2=kk+1;
                        end

                        if(hist(kk)>hist(k1)&&hist(kk)>hist(k2)...
                            &&hist(kk)>mag_thr)
                            bin=kk-1+0.5*(hist(k1)-hist(k2))/(hist(k1)+hist(k2)-2*hist(kk));
                            if(bin<0)
                                bin=HIST_BIN+bin;
                            elseif(bin>=HIST_BIN)
                                bin=bin-HIST_BIN;
                            end
                            key_number=key_number+1;
                            key_point_array(key_number,1)=k;
                            key_point_array(key_number,2)=j;
                            key_point_array(key_number,3)=i;
                            key_point_array(key_number,4)=(360/HIST_BIN)*bin;
                        end
                    end
                end
            end
        end
    end

    uni1=key_point_array(:,[1,2]);
    [~,ii,~]=unique(uni1,'rows','first');
    key_point_array=key_point_array(sort(ii)',:);

    constraint_number=size(key_point_array,1);
    if(key_point_array(constraint_number,1)==0&&key_point_array(constraint_number,2)==0)
        key_point_array=key_point_array(1:constraint_number-1,:);
    else
        key_point_array=key_point_array(:,:);
    end
end

function[hist,max_value]=Hist_Oritation(x,y,scale,gradient,angle,n)









    radius=round(6*scale);

    [MM,NN,~]=size(gradient);
    radius_x_left=x-radius;
    radius_x_right=x+radius;
    radius_y_up=y-radius;
    radius_y_down=y+radius;


    if(radius_x_left<=0)
        radius_x_left=1;
    end
    if(radius_x_right>NN)
        radius_x_right=NN;
    end
    if(radius_y_up<=0)
        radius_y_up=1;
    end
    if(radius_y_down>MM)
        radius_y_down=MM;
    end

    center_x=x-radius_x_left+1;
    center_y=y-radius_y_up+1;


    sub_gradient=gradient(radius_y_up:radius_y_down,radius_x_left:radius_x_right);
    sub_angle=angle(radius_y_up:radius_y_down,radius_x_left:radius_x_right);
    W=sub_gradient;
    bin=round(sub_angle*n/360);


    bin(bin>=n)=bin(bin>=n)-n;
    bin(bin<0)=bin(bin<0)+n;


    temp_hist=zeros(1,n);
    [row,col]=size(sub_angle);
    for i=1:1:row
        for j=1:1:col

            if(((i-center_y)^2+(j-center_x)^2)<=radius^2)
                temp_hist(bin(i,j)+1)=temp_hist(bin(i,j)+1)+W(i,j);
            end
        end
    end


    hist=zeros(1,n);
    hist(1)=(temp_hist(35)+temp_hist(3))/16+...
    4*(temp_hist(36)+temp_hist(2))/16+temp_hist(1)*6/16;
    hist(2)=(temp_hist(36)+temp_hist(4))/16+...
    4*(temp_hist(1)+temp_hist(3))/16+temp_hist(2)*6/16;


    hist(3:n-2)=(temp_hist(1:n-4)+temp_hist(5:n))/16+...
    4*(temp_hist(2:n-3)+temp_hist(4:n-1))/16+temp_hist(3:n-2)*6/16;

    hist(n-1)=(temp_hist(n-3)+temp_hist(1))/16+...
    4*(temp_hist(n-2)+temp_hist(n))/16+temp_hist(n-1)*6/16;
    hist(n)=(temp_hist(n-2)+temp_hist(2))/16+...
    4*(temp_hist(n-1)+temp_hist(1))/16+temp_hist(n)*6/16;


    max_value=max(hist);
end
