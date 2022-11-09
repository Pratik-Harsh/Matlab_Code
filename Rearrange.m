function [LD]=Rearrange(LD,ref)

i=find(LD(:,2)==ref);
if size(i,1)==0
    i=find(LD(:,3)==ref);
    for j=1:size(i,1)
        m=LD(i(j),2);
        LD(i(j),2)=LD(i(j),3);
        LD(i(j),3)=m;
    end
end
newLD=LD(i,:);
count=0;
LD(i,:)=[];

while size(LD,1)~=0    
    count=count+1;
    while count<=size(newLD,1)
        temp=find(LD(:,2)==newLD(count,3));
        if size(temp,1)~=0
            for k=1:size(temp,1)
                t=newLD;
                newLD=[t; LD(temp(k),:)];
            end
            LD(temp,:)=[];
        end
        count=count+1;
    end

    for i=1:size(LD,1)
        m=LD(i,2);
        LD(i,2)=LD(i,3);
        LD(i,3)=m;
    end
    i=find(LD(:,2)==ref);
    if size(i,1)~=0
        t=newLD;
        newLD=[t; LD(i,:)];
        LD(i,:)=[];
    end
    count=0;
end

[~,Y]=sort(newLD(:,1));

LD=newLD(Y,:);

