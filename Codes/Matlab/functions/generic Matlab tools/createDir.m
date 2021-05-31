function createDir(dirPath)

    %createDir(dirPath)
    %   creates a directory only if it doesn't exist. The reason is because
    %   mkdir works only if the directory doesn't exist. This function
    %   thus works for both the case where it exists or not
    %   
    %   INPUT:
    %       dirPath: path of the directory to create
    
    if ~exist(dirPath,'dir')
        mkdir(dirPath);
    end

end