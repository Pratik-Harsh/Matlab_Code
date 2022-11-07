clc;
clear all;

%% ---Input data of 33-bus test system--- %%

[BD,LD,TL]=data33(); %%---Bus, line, and tie data---%%
Vbase=12.66; %%---Base Voltage in kV---%%
Sbase=1; %%---Base Power in MVA---%%

%% ---Input data of 69-bus test system--- %%

[BD,LD,TL]=data69(); %%---Bus, line, and tie data---%%
Vbase=12.66; %%---Base Voltage in kV---%%
Sbase=1; %%---Base Power in MVA---%%

%% ---Input data of 84-bus test system--- %%

[BD,LD,TL]=data84(); %%---Bus, line, and tie data---%%
Vbase=11.40; %%---Base Voltage in kV---%%
Sbase=1; %%---Base Power in MVA---%%

%% ---Input data of 136-bus test system--- %%

[BD,LD,TL]=data136(); %%---Bus, line, and tie data---%%
Vbase=13.8; %%---Base Voltage in kV---%%
Sbase=1; %%---Base Power in MVA---%%

%% ---Input data of 417-bus test system--- %%

[BD,LD,TL]=data417(); %%---Bus, line, and tie data---%%
Vbase=10; %%---Base Voltage in kV---%%
Sbase=100; %%---Base Power in MVA---%%

%% ---Other informations--- %%

nbr=size(LD,1);
ntl=size(TL,1);
nbus=size(BD,1);

itermax=20;
iter=1;
Loss=zeros(itermax,1);

ref=find(BD(:,2)==3,1);

[Yb, ~] = Ybus(LD,nbr,nbus);
[Vmag, ~, Pcalc, ~]= NR_method(BD,Yb,nbus);

Vmin=min(Vmag);
Vmax=max(Vmag);

Loss(iter)=sum(Pcalc);
iter=iter+1;

x=sort(TL(:,1));
fprintf('Tie switch in initial Configuration: %s \n',num2str(x'))
fprintf('The total active power loss in the initial configuration of %d bus system is %f kW. \n\n', nbus, Loss(1)*Sbase*1000)

tstart=tic;

%% ---List of switches--- %%

[LD]=Rearrange(LD,ref);
[bibc, ~]=BIBC(LD,nbus,ref);

n_1=3;
n_2=2;

sw=zeros(nbr+ntl,1);
sw(TL(:,1))=sw(TL(:,1))+ones(ntl,1);

%%---Exclusion of Type-3 switch---%%
for i=1:ntl
    sw(LD(:,1))=sw(LD(:,1))+abs(bibc(:,TL(i,2))-bibc(:,TL(i,3)));
end
tempsw=sw;

%% ---Sequential Switch Opening using minimum current method--- %%

LD=[LD;TL];

t=0;

[Yb, A] = Ybus(LD,size(LD,1),nbus);
[Vmag, theta, Pcalc, ~]= NR_method(BD,Yb,nbus);

Vbus=Vmag.*(cos(theta)+1i.*sin(theta));
Vbr=Vbus(LD(:,2))-Vbus(LD(:,3));
Z=LD(:,4)+1i.*LD(:,5);
Ibr=Vbr./Z;

% disp(sum((abs(Ibr).^2).*LD(:,4)).*Sbase.*1000)

[~,index1]=sort((abs(Ibr)),'ascend');
index1=LD(index1,1);

t1=1;

while t<ntl
    if sw(index1(t1))>0
        index=find(LD(:,1)==index1(t1));
        tempA=A;
        tempA(index,:)=[];

        if rank(tempA)==nbus-1
            t=t+1;
            TL(t,:)=LD(index,:);
            LD(index,:)=[];

            A=tempA;
            t1=t1+1;
        else
            t1=t1+1;
        end
    else
        t1=t1+1;
    end

end

%% ---Branch exchange method using BIBC matrix--- %%

nbr=size(LD,1);
[LD]=Rearrange(LD,ref);
[bibc, bcbv]=BIBC(LD,nbus,ref);

%%---Exclusion of Type-1 switch---%%
temp=sum(bibc,1);
for i=1:nbus
    if temp(i)<n_1
        for j=1:nbr
            if bibc(j,i)==1
                sw(LD(j,1))=0;
            end
        end
    end
end

%%---Exclusion of Type-2 switch---%%
for i=1:nbr
    x=find(LD(:,2)==LD(i,3));
    if size(x,1)>n_2
        sw(LD(i,1))=0;
    end
end

Vbus=BD(:,3).*(cos(BD(:,4))+1i.*sin(BD(:,4)));

Sinj=(BD(:,5)-BD(:,7))+1i*(BD(:,6)-BD(:,8));

[Vbus,Iinj]= dlf(ref,nbus,Sinj,Vbus,bibc,bcbv);
Ibr=-bibc*Iinj;
Vmag=abs(Vbus);

Loss(iter)=sum((abs(Ibr).^2).*LD(:,4));

iter=iter+1;
maxdPloss=inf;
iter1=0;

while iter1<1

    %%---Selection of switch pair for the exchange operation---%%

    maxdPloss=-inf;

    for k1=1:ntl

        m1=TL(k1,2);
        n1=TL(k1,3);
        Rk1=TL(k1,4);

        loop=abs(bibc(:,m1)-bibc(:,n1));

        for k2=1:nbr
            if loop(k2)==1 && sw(LD(k2,1))>0

                Ibrnew=Ibr+loop.*(1-2.*(bibc(:,m1).*bibc(k2,m1)+bibc(:,n1).*bibc(k2,n1))).*Ibr(k2);

                dPloss=-sum((abs(Ibrnew).^2).*LD(:,4))+sum((abs(Ibr).^2).*LD(:,4))-abs(Ibr(k2)^2)*Rk1;

                if maxdPloss<dPloss
                    maxdPloss=dPloss;
                    index2=k2;
                    index1=k1;
                end
            end

        end
    end

    %%---Updating the configuration---%%

    tempLD=[LD; TL(index1,:)];
    temp=tempLD(index2,:);
    tempLD(index2,:)=[];

    [tempLD]=Rearrange(tempLD,ref);
    [tempbibc, tempbcbv]=BIBC(tempLD,nbus,ref);

    [tempVbus,Iinj]= dlf(ref,nbus,Sinj,Vbus,tempbibc,tempbcbv);

    tempIbr=-tempbibc*Iinj;
    temploss=sum((abs(tempIbr).^2).*tempLD(:,4));
    tempVmag=abs(Vbus);

    if temploss<Loss(iter-1) && min(tempVmag)>=Vmin && max(tempVmag)<=Vmax

        LD=tempLD;

        TL(index1,:)=[];
        TL(ntl,:)=temp;

        Vbus=tempVbus;
        bibc=tempbibc;
        Ibr=tempIbr;

        sw=tempsw;

        %%---Exclusion of Type-1 switch---%%
        temp=sum(bibc,1);
        for i=1:nbus
            if temp(i)<n_1
                for j=1:nbr
                    if bibc(j,i)==1
                        sw(LD(j,1))=0;
                    end
                end
            end
        end

        %%---Exclusion of Type-2 switch---%%
        for i=1:nbr
            x=find(LD(:,2)==LD(i,3));
            if size(x,1)>n_2
                sw(LD(i,1))=0;
            end
        end

        Loss(iter)=temploss;
        iter=iter+1;

    else
        Loss(iter)=Loss(iter-1);
        iter=iter+1;

        temp=TL(index1,:);
        TL(index1,:)=[];
        TL(ntl,:)=temp;
        iter1=iter1+1;

    end

end

Loss(iter:end)=[];
Loss=Loss.*Sbase*1000;

tstop=toc(tstart);

%% ---RESULTS--- %%

x=sort(TL(:,1));
fprintf('Tie switch in final Configuration: %s \n',num2str(x'))
fprintf('The total active power loss in the final configuration of %d bus system is %f kW. \n\n', nbus, Loss(end))

fprintf('The total computational time is %f sec. \n\n',tstop)









