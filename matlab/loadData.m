
function [data,messages] = loadData(path)
%function data = loadData(path)
%data.x = [x1,x2,x3,x4,...]
%data.y = [y1,y2,y3,y4,...]
%data.data = [d_11,d_12,d_13,d_14;d_21,d_22,d_23,d24;...]
%Works with binary or delimited ascii/unicode data. Save the data file as a
%.txt file if you want to read it in as ascii or unicode. Otherwise, it
%will default to binary.

try
    fidP = fopen(path,'r','ieee-le.l64');
catch
    data = [];
    messages = {'DATA FILE FAILED TO LOAD'};
    return
end
if fidP == -1
    data = [];
    messages = {'BAD FILE PATH'};
    return
end

try
    index = 0;
    fseek(fidP,0,'eof');
    fileEnd = ftell(fidP);
    fseek(fidP,0,'bof');
    while ftell(fidP) < fileEnd - 8
        index = index + 1;
        messages = [];
        if strcmp(path(end-3:end),'.txt')
            %ascii or unicode
            xNum(index) = fscanf(fidP,'%e',1);
            yNum(index) = fscanf(fidP,'%e',1);
            data(index).x = fscanf(fidP,'%e',xNum(index));
            data(index).y = fscanf(fidP,'%e',yNum(index));
            if yNum(index) == 0
                data(index).data = fscanf(fidP,'%e',[xNum(index),1]);
                if size(data.data,1)*size(data.data,2) ~= xNum(index)*yNum(index)
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'};...
                        {'It might not like how your decimal numbers are formatted'}];
                end
            else
                data(index).data = fscanf(fidP,'%e',[xNum(index),yNum(index)]);
                if size(data.data,1)*size(data.data,2) ~= xNum(index)*yNum(index)
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'};...
                        {'It might not like how your decimal numbers are formatted'}];
                end
            end
            if isempty(xNum)
                messages ={'DATA MISLOADED. xNum is empty...'};
            end
        else
            %binary
            xNum(index) = fread(fidP,1,'float64');
            if(xNum(1) < 1) %If it looks like the file is wrong endian,
                %close the file and re-open it with little endian.
                fclose(fidP);
                fidP = fopen(path,'r','ieee-be.l64');
                xNum(index) = fread(fidP,1,'float64');
            end
            yNum(index) = fread(fidP,1,'float64');
            data(index).x = fread(fidP,xNum(index),'float64');
            data(index).y = fread(fidP,yNum(index),'float64');
            if or(yNum(index) == 0, yNum(index) == 1)
                data(index).data = fread(fidP,[xNum(index),1],'float64');
                if size(data(index).data,1)*1 ~= xNum(index)*1
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'}];
                end
            else
                data(index).data = fread(fidP,[xNum(index),yNum(index)],'float64');
                if size(data(index).data,1)*size(data(index).data,2) ~= xNum(index)*yNum(index)
                    messages =[{'DATA MISLOADED (size != xNum*yNum)'};...
                        {'Is your file formatted correctly?'}];
                end
            end
            if isempty(xNum)
                messages =[{'DATA MISLOADED'};...
                    {'Is your file opening a big endian (PowerPC) file'};...
                    {'file on a little endian (x86) machine?'}];
            end
        end
    end
    
    fclose(fidP);
catch
    data = [];
    messages = {'DATA FAILED TO LOAD'};
    fclose(fidP);
    return
end
