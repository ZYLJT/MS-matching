function P2fine = Finepoints(P1,P2)%��ƥ���������Խ���"subpixel"����ľ�ϸƥ��

    len=length(P1);
    affmat=fitgeotrans(P1,P2,'projective');
    P2pro=[P1,ones(len,1)]*affmat.T;
    P2pro(:,1)=P2pro(:,1)./P2pro(:,3);
    P2pro(:,2)=P2pro(:,2)./P2pro(:,3);
    P2pro(:,3)=[];
    devia_P=(P2-P2pro).^2;
    devia_P=sqrt(sum(devia_P,2));%��������P2��P2pro֮���ƫ��devia_P��������ת��Ϊŷ����þ��롣
    max_Devia=max(devia_P);
    iteration=0;
    P2fine=P2;
    while max_Devia>0.05&&iteration<10%��ֹ����
        iteration=iteration+1;
        [~,index]=sort(devia_P);
        %ind1=round(1/5*length(index));
        ind1=round(1/4*length(index));
        P2fine(index(ind1:end),:)=P2pro(index(ind1:end),:);
        affmat=fitgeotrans(P1,P2fine,'projective');
        P2pro=[P1,ones(len,1)]*affmat.T;
        P2pro(:,1)=P2pro(:,1)./P2pro(:,3);
        P2pro(:,2)=P2pro(:,2)./P2pro(:,3);
        P2pro(:,3)=[];
        devia_P=(P2fine-P2pro).^2;
        devia_P=sqrt(sum(devia_P,2));
        max_Devia=max(devia_P);
    end
end
