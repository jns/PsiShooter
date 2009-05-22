function writeFile(X,Y,Data,path, units)
%writeFile(X,Y,Data,fileName)

if Y == 0
    yLength = 1;
else
    yLength=length(Y);
end

% X and Y must be in cm. 
% data can be in ev or ergs to convert output to ergs
conv = 1;
if (strcmp(units, 'ev')) 
    conv = 1.60217646e-12;
end

try
    fidP = fopen(path,'w','ieee-le.l64');
    if strcmp(path(end-3:end),'.txt')
        %ascii or unicode
        fprintf(fidP,'%e\n',length(X));
        fprintf(fidP,'%e\n',yLength);
        fprintf(fidP,'%e ',X);
        %fprintf(fidP,'%s','\n');
        fprintf(fidP,'%e ',Y);
        %fprintf(fidP,'%s','\n');
        fprintf(fidP,'%e ',Data*conv);
    else
        %binary
        fwrite(fidP,length(X),'float64');
        fwrite(fidP,yLength,'float64');
        fwrite(fidP,X,'float64');
        fwrite(fidP,Y,'float64');
        fwrite(fidP,Data*conv,'float64');
    end
    
    fclose(fidP);
catch
    fclose(fidP);
    return
end