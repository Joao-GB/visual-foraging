function filename = fileVersion(baseName, extension)
    filename = [baseName, extension];
    counter = 1;
    
    while exist(filename, 'file')
        filename = sprintf('%s(%d)%s', baseName, counter, extension);
        counter = counter + 1;
    end
end