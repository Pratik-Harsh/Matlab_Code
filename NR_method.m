function [Vmag, theta, Pcalc, Qcalc]= NR_method(BD,Yb,nbus)

Vmag=BD(:,3); %%---Bus Voltage Magnitude
theta=BD(:,4); %%---Bus Voltage Angle
Ymag=abs(Yb); 
Yang=angle(Yb);

Pcomp=zeros(nbus,1); %%----Active power of the composite load
Qcomp=zeros(nbus,1); %%----Reactive power of the composite load
c=[1 0 0];
d=[1 0 0];

npv=0;
for i=1:nbus
    if BD(i,2)==2
        npv=npv+1; %%---Number of PV bus
    end
end

iter=0;
itermax=200;
mis=ones(2*(nbus-1),1);
tol=10^(-8);

while max(abs(mis))>tol
    iter=iter+1;

    Pcalc=zeros(nbus,1);
    Qcalc=zeros(nbus,1);
    J=zeros(2*(nbus-1),2*(nbus-1)); %%---Initialisation of Jacobian Matrix
    m=0;
    p=nbus-1+m;
    for i=1:nbus
        for j=1:nbus
            Pcalc(i)=Pcalc(i)+Vmag(i)*Vmag(j)*Ymag(i,j)*cos(theta(i)-theta(j)-Yang(i,j)); %%---Calculated value of Bus Active Power Injection
            Qcalc(i)=Qcalc(i)+Vmag(i)*Vmag(j)*Ymag(i,j)*sin(theta(i)-theta(j)-Yang(i,j)); %%---Calculated value of Bus Reactive Power Injection
        end
        Pcomp(i)=BD(i,7)*(c(1)+c(2)*Vmag(i)+c(3)*Vmag(i)*Vmag(i));
        Qcomp(i)=BD(i,8)*(d(1)+d(2)*Vmag(i)+d(3)*Vmag(i)*Vmag(i));

        if BD(i,2)~=3
            m=m+1;
            n=0;
            q=nbus-1+n;
            for j=1:nbus
                if BD(j,2)~=3
                    n=n+1;
                    if i==j
                        J(m,n)=-Qcalc(i)-Vmag(i)*Vmag(i)*imag(Yb(i,i));
                    else
                        J(m,n)=Vmag(i)*Vmag(j)*Ymag(i,j)*sin(theta(i)-theta(j)-Yang(i,j));
                    end
                    if BD(j,2)==1
                        q=q+1;
                        if i==j
                            J(m,q)=Pcalc(i)+Vmag(i)*Vmag(i)*real(Yb(i,j));
                        else
                            J(m,q)=Vmag(i)*Vmag(j)*Ymag(i,j)*cos(theta(i)-theta(j)-Yang(i,j));
                        end
                    end
                end
            end
            
            if BD(i,2)==1
                p=p+1;
                n=0;
                q=nbus-1+n;
                for j=1:nbus
                    if BD(j,2)~=3
                        n=n+1;
                        if i==j
                            J(p,n)=Pcalc(i)-Vmag(i)*Vmag(i)*real(Yb(i,i));
                        else
                            J(p,n)=-Vmag(i)*Vmag(j)*Ymag(i,j)*cos(theta(i)-theta(j)-Yang(i,j));
                        end
                        if BD(j,2)==1
                            q=q+1;
                            if i==j
                                J(p,q)=Qcalc(i)-Vmag(i)*Vmag(i)*imag(Yb(i,j));
                            else
                                J(p,q)=Vmag(i)*Vmag(j)*Ymag(i,j)*sin(theta(i)-theta(j)-Yang(i,j));
                            end
                        end
                    end
                end
            end
        end
    end
    
    Psp=BD(:,5)-Pcomp; %%---Specified value of Bus Active Power Injection
    Qsp=BD(:,6)-Qcomp; %%---Specified value of Bus Reactive Power Injection
    
    m=0;
    n=0;
    misP=zeros(nbus-1,1);
    misQ=zeros(nbus-1,1);
    for i=1:nbus
        if BD(i,2)~=3
            m=m+1;
            misP(m,1)=Psp(i,1)-Pcalc(i,1);
            if BD(i,2)==1
                n=n+1;
                misQ(n,1)=Qsp(i,1)-Qcalc(i,1);
            end
        end
    end
    
    mis=[misP; misQ]; %%---Mismatch Vector    
        
    del=J\mis;

    m=0;
    n=nbus-1;
    for i=1:nbus
        if BD(i,2)~=3
            m=m+1;
            theta(i)=theta(i)+del(m); %%---Updated Value of Bus Voltage Angle
            if BD(i,2)==1
                n=n+1;
                Vmag(i)=Vmag(i)+Vmag(i)*del(n); %%---Updated Value of Bus Voltage Magnitude
            end
        end
    end
    if iter==itermax
        disp('Load Flow did not converged');
        break;
    end
end
