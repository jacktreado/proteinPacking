%% function to make easy .txt files out of hydrogenated PDB file

function cleanSinglePDB(pdbstr,radLoc,outputstr)

% load atom sizes info
if radLoc(end) == '/'
    radLoc(end) = [];
end
radFid = fopen([radLoc '/radii.dat']);

% get radii data
radData         = textscan(radFid,'%s %s %s %f');
radResName      = radData{1};
radAtomName     = radData{2};
radBondType     = radData{3};
radValue        = radData{4};

% load pdb file
pdbStruct       = pdbread(pdbstr);

% get atomic coordinates
pdbAtoms        = pdbStruct.Model.Atom;

resNames        = {pdbAtoms.resName}';
resSeq          = cell2mat({pdbAtoms.resSeq}');
atomNames       = {pdbAtoms.AtomName}';
chainNames      = {pdbAtoms.chainID}';
ax              = cell2mat({pdbAtoms.X}');
ay              = cell2mat({pdbAtoms.Y}');
az              = cell2mat({pdbAtoms.Z}');
bf              = cell2mat({pdbAtoms.tempFactor}');

% loop over atoms, assign radii, print to output file string
outFid = fopen(outputstr,'w');

% output header
NA = length(ax);
N = max(resSeq);
fprintf(outFid,'%d\n',NA);
fprintf(outFid,'%d\n',N);
for aa = 1:NA
    % get residue index for radius lookup
    resNameIndex    = strcmp(resNames{aa},radResName);
    atomNameIndex   = strcmp(atomNames{aa},radAtomName);
    radiusIndex     = resNameIndex & atomNameIndex;
    
    % get radius value
    radtmp          = radValue(radiusIndex);
    if sum(radiusIndex) > 1     
        radtmp      = radtmp(1);
    elseif sum(radiusIndex) == 0
        atomFirst = atomNames{aa};
        atomFirst = atomFirst(1);
        switch atomFirst
            case 'H'
                radtmp = 1.1;
            case 'O'
                radtmp = 1.4;
            case 'N'
                radtmp = 1.3;
            case 'C'
                radtmp = 1.5;
            case 'S'
                radtmp = 1.75;
            otherwise
                radtmp = 1.4;
        end
    end
    
    % print info to file
    fprintf(outFid,'%6d %10s %10s %10s %10f %20.8f %20.8f %20.8f\n',...
        resSeq(aa),resNames{aa},atomNames{aa},chainNames{aa},radtmp,ax(aa),ay(aa),az(aa));
%     fprintf('%6d %10s %10s %10s %10f %20.8f %20.8f %20.8f\n',...
%         resSeq(aa),resNames{aa},atomNames{aa},chainNames{aa},radtmp,ax(aa),ay(aa),az(aa));
end




% close files
fclose(outFid);
fclose(radFid);
end

