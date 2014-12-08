function [Ar ] = mylatt2(file1,file2,newfile)
P=dlmread(file1)';
Q=dlmread(file2)';
    for i=1:4
       for j=1:4
                Cp(i,j)=floor(P(i,j))==P(i,j);
                Cq(i,j)=floor(Q(i,j))==Q(i,j);
       end
    end
    if isequal(Cp,ones(4))==0
        disp('No integer entries')
        Ar=0;
    elseif isequal(Cq,ones(4))==0
          disp('No integer entries')
          Ar=0;
    elseif isequal(rank(P),4)==0
          disp('The rank is not 4') 
          Ar=0;
    elseif isequal(P(1,:),ones(1,4))==0
        disp('Non-homogeneous coordinates')
        Ar=0;
    elseif isequal(Q(1,:),ones(1,4))==0
        disp('Non-homogeneous coordinates')
        Ar=0;
  else 
    v=[1;2;3;4];
    Per2=perms(v);
    q=1;
    for i=1:24
        a=Per2(i,1);
        b=Per2(i,2);
        c=Per2(i,3);
        d=Per2(i,4);
        y=[P(:,a),P(:,b),P(:,c),P(:,d)];
        if y==Q
            I=eye(4);
            Z1=[I(:,a),I(:,b),I(:,c),I(:,d)];
            Ar{q}=Z1;
            q=q+1;
        end
        Anew=Q*inv(y);
            for i=1:4
                for j=1:4
                aij=Anew(i,j);
                r= str2num(sprintf('%.4f',aij));
                T(i,j)=floor(r)==r;
                end
            end   
            if isequal(T,ones(4))==1
            Ar{q}=round(Anew); %because of some computations sometimes there are numbers as -5.5511e-17
            q=q+1;
            end          
    end
    le=size(Ar);
    dlmwrite(newfile,Ar{1},'delimiter',' ');
    for i=2:le(2)
        dlmwrite(newfile,'-----','-append','Delimiter','');
        dlmwrite(newfile,Ar{i},'-append','Delimiter',' ');
    end
    end
        
    return
    
%   It returns a set of matrices encoding an affine map that sends P to Q in a new file. 


end

