function [ c ] = orbitsize( X,n,v )
tic
if nargin <1 
    error('orbitsize :  X is a required input')
end
if nargin < 2
 error('orbitsize :  n is a required input')
end
%Determine the generator:
if X=='A'
    An=[];
    for i=1:n
        for j=1:n+2
            if j-i==1
                An(i,j)=1;
            elseif j-i==2
                An(i,j)=-1;
            else An(i,j)=0;
            end
        end
    end
    G=An;
end

if X=='D'
    Dn=[];
    for i=1:n-1
        for j=1:n+1
            if j-i==1
                Dn(i,j)=1;
            elseif j-i==2
                Dn(i,j)=-1;
            else Dn(i,j)=0;
            end
        end
    end
    for j=1:(n-1)
            Dn(n,j)=0;
    end
    Dn(n,n)=1;
    Dn(n,n+1)=1;
G=Dn;
end
 if X=='E'
     if isequal(n,8)==0
         error('orbit; for generator E dimension must be 8')
     else E8=[];
         for i=1:6
             for j=1:9
                 if j-i==1
                E8(i,j)=1;
            elseif j-i==2
                E8(i,j)=-1;
                 end
             end
         end
          for j=1:6
         E8(7,j)=0;
          end
          E8(7,7)=1;
          E8(7,8)=1;
          E8(8,8)=0;
          E8(8,1)=0;
          for j=2:9
              E8(8,j)=1/2;
          end
      G=E8;    
     end
 end
 
if X=='F'
     if isequal(n,4)==0
         error('orbit; for generator F n must be 4')
     else F4=[];
         for i=1:2
             for j=1:4
                 if j-i==1
                F4(i,j)=1;
            elseif j-i==2
                F4(i,j)=-1;
                 end
             end
         end
          for j=[1,2,3,5]
             F4(3,j)=0;
          end
          F4(3,4)=1;
          for j=2:5
          F4(4,j)=-1/2;
          end
          G=F4;
     end
 end
          
if X=='H'
        if isequal(n,3)==0 & isequal(n,4)==0
         error('orbit; for generator H n must be 3 or 4')
        elseif n==3 
         tau=1/2+1/2*sqrt(5);
         H3=[0,2,0,0;0,-(tau+1),1,-tau;0,0,0,2]; 
         G=H3;
        elseif n==4
            tau=1/2+1/2*sqrt(5);
            H4=[0,-(tau+1),1,1,1;0,-1,1,0,0;0,0,-1,1,0;0,0,0,-1,1];
            G=H4;
        end
         if nargin < 3
             v(1)=0;
                v(2)=1;
                 for i=3:n+1
                v(i)=i-1;
                 end
         end
end

 
%define a vector v out of evary hyperplane if you dont recive the input v.
lg=size(G);
if nargin < 3 & X~='H'
    v(1)=0;
    v(2)=1;
    for i=3:lg(2)
        v(i)=1;
    end
        t=1;
        mu=G*v';
    while or(any(mu==0),length(unique(mu))~=length(mu))
        t=t+1;
        for i=1:lg(2)-2
        v(i+2)=t.^i;
        end
        mu=G*v';
    end
end


%%%calculates all the possible images%%
R(1,:)=v;
le=size(R);
i=1;
while i<=le(1)
    for j=1:lg(1)
        a=G(j,:);
        r=R(i,:);
        reff=r-2*(((r*a')/(a*a'))*a);
        ref=round(reff*10000)/10000;
        if all(ismember(R,ref,'rows')==0)
            R=[R;ref];
            le=size(R);
        end
    end
    i=i+1;
end
c=le(1);
toc
end



%NOTES: As the program has problems comparing vectors with long decimals, when a vector v is not given a defined it appart for the group H since computing a vector in the way I did it for the others was not working right becuse of the comparisong of decimals.



