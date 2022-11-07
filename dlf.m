function [Vbus,Iinj]= dlf(ref,nbus,Sinj,Vbus,bibc,bcbv)

delV=1;
Vmag=abs(Vbus);
tol=10e-8;

while delV>tol
    Iinj=conj(Sinj./Vbus);

    V=Vbus(ref)+bcbv*bibc*Iinj;

    for i=1:nbus
        if i<ref
            Vbus(i)=V(i);
        elseif i>ref
            Vbus(i)=V(i-1);
        end
    end

    delV=max(abs(Vmag-abs(Vbus)));
    Vmag=abs(Vbus);
end
