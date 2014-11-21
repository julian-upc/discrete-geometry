function [ t] = test_chi( file1,file2,result)
Ar=mylatt2(file1,file2,'myresult.txt');
R=dlmread(result)';
le=size(Ar);
t=1;
i=1;
    while isequal(Ar{i},R)==0;
        i=i+1;
        if i>le(2)
            t=0;
            break
        else continue
        end  
    end  
%The functions returns one if there the matrix in the file result is an
%afffine map sending one matrix into another and 0 otherwise.


end

