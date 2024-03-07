function[y0]=LIP(xvec,yvec,x0);

    
    n=length(xvec);
    m=length(x0);
    y0=zeros(m,1);

    for k=1:m;
        if x0(k)==xvec(1);
            y0(k)=yvec(1); 
        elseif x0(k)==xvec(n);
            y0(k)=yvec(n);
        else;
            z = zeros(length(xvec),1);
            z(xvec<=x0(k))=1;
            j=sum(z); % this determines the lower bracket for x0 @
        y0(k)=yvec(j)+((yvec(j+1)-yvec(j))/(xvec(j+1)-xvec(j)))*(x0(k)-xvec(j));
        end;
    end;
    
