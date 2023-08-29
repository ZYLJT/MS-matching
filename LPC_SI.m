function[or,Gradient]=LPC_SI(im,C,scales,norient,B)

    if(~exist('C'))
        C=2;
    end
    if(~exist('scales'))
        scales=[1,3/2,2];w=[1,-3,2];
    end
    if(~exist('norient'))
        norient=6;
    end
    nscale=length(scales);
    [row,col]=size(im);

    if(~exist('B'))
        B=round(min(row,col)/16);
    end
    im=double(im);
    [filter,Orientation]=LPSO_logGabor_2D(im,norient,nscale,scales,0.55);
    imfft=fft2(im);
    s_lpc=ones(row,col,norient);

    for o=1:norient
        for s=1:nscale
            M(:,:,s)=ifft2(imfft.*filter{s,o});
            s_lpc(:,:,o)=s_lpc(:,:,o).*(M(:,:,s).^w(s));
        end
        e=abs(M(:,:,1));
        e_center=e(B+1:end-B,B+1:end-B);
        e_mean=mean(e_center(:));
        e_std=std(e_center(:));
        T=e_mean+2*e_std;%ÔëÉùãÐÖµ
        %T=log(1+e_mean/e_std);%ÔëÉùãÐÖµ
        %T = log(1+e_std^2/e_mean^2);
        %e=max(0,e-T);
        [row,col]=size(e);
        for i=1:row
            for j=1:col
                if e - T>=0
                    e(i,j) = e(i,j)- T*exp(T-e(i,j)); 
                else
                     e(i,j)=0;
                end
            end
        end
        
        
        energy(:,:,o)=e;

    end
    s_lpc_map=cos(angle(s_lpc));
    s_lpc_map(s_lpc_map<0)=0;
    lpc_map=(sum(s_lpc_map.*energy,3))./(sum(energy,3)+C);
%     for j=1:o%6
%         s_lpc_map=cos(angle(s_lpc));
%         s_lpc_map(s_lpc_map<0)=0;
%         lpc_map=(sum(s_lpc_map.*energy,3))./(sum(energy,3)+C);
%     end

    or=Orientation;
    Gradient=lpc_map;
    return;
