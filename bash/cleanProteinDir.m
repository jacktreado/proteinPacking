function cleanProteinDir(inDir,outDir,radLoc,frmt)
%% FUNCTION to process pdb files with added hydrogens

% remove / at end of inDir str
if inDir(end) == '/'
    inDir(end) = [];
end

if outDir(end) == '/'
    outDir(end) = [];
end

% get list of files
flist = dir([inDir '/*' frmt]);
NF = length(flist);
if NF == 0
    fprintf('No files found in directory %s with format %s\n',inDir,frmt);
    error('no files found with input format, killing program.');
end

% loop over files, clean pdb and print to new file
fprintf('Looping over NF = %d files...\n',NF);
for ff = 1:NF        
    % get current file
    pdbname = flist(ff).name;
    pdbdir = flist(ff).folder;
    pdbcode = pdbname(1:end-4);
    pdbf = [pdbdir '/' pdbname];
    
    % print start of file processing
    fprintf(' On file ff = %d/%d, pdb file = %s\n',ff,NF,pdbname);
    
    % get name of output file
    outputf = [outDir '/' pdbcode '.dat'];
    
    % clean single pdb
    % note, assuming radii file is located in working directory
    fprintf('-- Printing data in %s to %s...',pdbname,outputf);
    cleanSinglePDB(pdbf,radLoc,outputf);
    fprintf('.Done!\n');
end

end