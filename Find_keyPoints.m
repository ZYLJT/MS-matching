function key_point_array=Find_keyPoints...
    (Corner_function,...
    sigma,...
    ratio,...
    first_points,...
    gradient,...
    angle,...
    Scale_Invariance)

    if nargin<8
        nOctaves=3;
    end


    HIST_BIN=36;
    SIFT_ORI_PEAK_RATIO=0.9;
    key_number=0;
    Fig=0;
    if(strcmp(Scale_Invariance,'YES'))
        Layer_points=4500;
        key_num=Layer_points;
        layers=1;
    else
        Layer_points=first_points;
        key_num=nOctaves*Layer_points;
        layers=nOctaves;
    end


    key_point_array=zeros(key_num,6);
    for i=1:1:layers
        Corner_temp_current=Corner_function{i};
        gradient_current=gradient{i};
        angle_current=angle{i};


        Fig=i;

        if Fig==1
            Diff='region';
        elseif Fig==2
            Diff='edge';
        else
            Diff='sharpedge';
        end
        a=max(Corner_temp_current(:));b=min(Corner_temp_current(:));Corner_temp_current=(Corner_temp_current-b)/(a-b);
        %corners_points=detectKAZEFeatures(Corner_temp_current,'Diffusion',Diff);
        corners_points = detectFASTFeatures(Corner_temp_current,'MinContrast',0.05);
        corners=corners_points.selectStrongest(Layer_points);
        NUMcorners=corners.Count;
        Pointscorners=corners.Location;
        Pointsstrength=corners.Metric;

        for L=1:NUMcorners
            k=Pointscorners(L,1);
            j=Pointscorners(L,2);
            strength=Pointsstrength(L,1);
            scale=sigma*ratio^(i-1);

            [hist,max_value]=Hist_Oritation(k,j,scale,gradient_current,angle_current,HIST_BIN);
            mag_thr=max_value*SIFT_ORI_PEAK_RATIO;
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
                    key_point_array(key_number,4)=((360/HIST_BIN)*bin)/2;
                    key_point_array(key_number,5)=strength;
                    key_point_array(key_number,6)=scale/sigma;
                end
            end
        end
    end

    uni1=key_point_array(:,[1,2]);
    [~,i,~]=unique(uni1,'rows','first');
    key_point_array=key_point_array(sort(i)',:);
    constraint_number=size(key_point_array,1);
    if(key_point_array(constraint_number,1)==0&&key_point_array(constraint_number,2)==0)
        key_point_array=key_point_array(1:constraint_number-1,:);
    else
        key_point_array=key_point_array(:,:);
    end

