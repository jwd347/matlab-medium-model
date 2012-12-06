fid=fopen('thermo.inp');
clc
ctLine=0;
while ctLine<2500
    
    ctLine=ctLine+1;
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    %     disp(tline)
    if tline(1)=='!'
        continue
    end
    if strfind(tline,'thermo')
        tline = fgetl(fid);
        continue
    end
    if strfind(tline,'END')
        %         tline = fgetl(fid);
        continue
    end
    %% Must have got to a real section
    str.name =  sscanf(tline(1:16),'%s');
    str.FullName = str.name;
    if str.name(end)=='+'
        str.name=[str.name(1:end-1) 'plus'];
    elseif str.name(end)=='-'
        str.name=[str.name(1:end-1) 'minus'];
    else
        
    end
    ind=regexp(str.name,'[()]');
    str.name(ind)='b';
    ind=regexp(str.name,'\W');
    str.name(ind)='_';
    if regexp(str.name(1),'[0-9]')
        str.name=['num_' str.name];
    end
    
    str.comments = tline(19:end);
    tline = fgetl(fid);
    str.ctTInt = str2num(tline(1:2));
    str.txRefCode = tline(4:9);
    str.txForumla = tline(11:50);
    str.swtCondensed = str2num(tline(51:52));
    str.mm = str2num(tline(53:65));
    str.Hf0 = str2num(tline(66:80));
    
    if str.ctTInt ==0
        
        tline = fgetl(fid);
        str.tRange = str2num(tline(1:22));
        str.Hf298Del0 = str2num(tline(66:end));
%         warning(' here')
    end
    for ctInterval=1:str.ctTInt
        tline = fgetl(fid);
        str.tRange{ctInterval} = str2num(tline(1:22));
        str.Hf298Del0{ctInterval} = str2num(tline(66:end));
        
        tline = fgetl(fid);
        a1 = str2num(tline(1:16));
        a2 = str2num(tline((1:16)+16));
        a3 = str2num(tline((1:16)+32));
        a4 = str2num(tline((1:16)+48));
        a5 = str2num(tline((1:16)+64));
        tline = fgetl(fid);
        a6 = str2num(tline(1:16));
        a7 = str2num(tline((1:16)+16));
        %a8 =blank
        b1 = str2num(tline((1:16)+48));
        b2 = str2num(tline((1:16)+64));
        str.a{ctInterval}=[a1 a2 a3 a4 a5 a6 a7];
        str.b{ctInterval}=[b1 b2];
    end
    
    
    %     disp(tline)
    disp( [str.name ' : ' str.FullName ' : ' str.comments] )
    strMaster.(str.name)=str;
    clear str
end

fclose(fid);
