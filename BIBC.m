function [bibc, bcbv]=BIBC(LD,nbus,ref)

nbr=size(LD,1);

bibc=zeros(nbr,nbus);
bcbv=zeros(nbus,nbr);

for i=1:nbr

    bibc(i,LD(i,3))=1;
    bcbv(LD(i,3),i)=LD(i,4)+1i*LD(i,5);

    br=LD(i,3);
    count=1;

    while count<=size(br,1)
        n=br(count,1);
        temp=find(LD(:,2)==n);
        if size(temp,1)~=0
            for j=1:size(temp,1)
                t=br;
                br=[t; LD(temp(j),3)];
                bibc(i,LD(temp(j),3))=1;
                bcbv(LD(temp(j),3),i)=LD(i,4)+1i*LD(i,5);
            end
        end
        count=count+1;
    end
end

bcbv(ref,:)=[];
