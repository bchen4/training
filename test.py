import subprocess

#specifying a list of parameters to run bowtie with
#each tuple contains
# *the mismatch threshold below which a site is considered a potential off-target (higher is more stringent)
# *the number of sites allowed (1 is minimum since each sgRNA should have one true site in genome)
# *the genome index against which to align the sgRNA sequences; these can be custom built to only consider sites near TSSs
# *a name for the bowtie run to create appropriately named output files
alignmentList = [(39,1,'large_data_files/hg19.ensemblTSSflank500b','39_nearTSS'),
                (31,1,'large_data_files/hg19.ensemblTSSflank500b','31_nearTSS'),
                (21,1,'large_data_files/hg19_maskChrMandPAR','21_genome'),
                (31,2,'large_data_files/hg19.ensemblTSSflank500b','31_2_nearTSS'),
                (31,3,'large_data_files/hg19.ensemblTSSflank500b','31_3_nearTSS')]

alignmentColumns = []
for btThreshold, mflag, bowtieIndex, runname in alignmentList:
    alignedFile = 'temp_bowtie_files/' + runname + '_aligned.txt'
    unalignedFile = 'temp_bowtie_files/' + runname + '_unaligned.fq'
    maxFile = 'temp_bowtie_files/' + runname + '_max.fq'
    bowtieString = '/Users/bbc/Tools/bowtie-1.1.2/bowtie -n 3 -l 15 -e '+str(btThreshold)+' -m ' + str(mflag) + ' --nomaqround -a --tryhard -p 16 --chunkmbs 256 ' + bowtieIndex + ' --suppress 5,6,7 --un ' + unalignedFile + ' --max ' + maxFile + ' '+ ' -q '+fqFile+' '+ alignedFile
    print bowtieString
    print subprocess.call(bowtieString, shell=True) #0 means finished without error
    #parse through the file of sgRNAs that exceeded "m", the maximum allowable alignments, and mark "True" any that are found
    try:
        with open(maxFile) as infile:
            sgsAligning = set()
            for i, line in enumerate(infile):
                if i%4 == 0: #id line
                    sgsAligning.add(line.strip()[1:])
    except IOError: #no sgRNAs exceeded m, so no maxFile created
        sgsAligning = set()                    
    alignmentColumns.append(libraryTable_new.apply(lambda row: row.name in sgsAligning, axis=1))
    
#collate results into a table, and flip the boolean values to yield the sgRNAs that passed filter as True
alignmentTable = pd.concat(alignmentColumns,axis=1, keys=zip(*alignmentList)[3]).ne(True)

#The big 3-layer loop
for (gene, transcript), group in v2Groups:
    geneSgIds = []
    geneLeftPositions = []
    empiricalSgIds = dict() 
    stringency = 0  
    while len(geneSgIds) < sgRNAsToPick and stringency < len(offTargetLevels):
        for sgId_v2, row in group.sort_values(('predicted score','predicted score'), ascending=False).iterrows():
            oligoSeq = upstreamConstant + row[('library table v2','sequence')] + downstreamConstant
            leftPos = row[('sgRNA info', 'position')] - (23 if row[('sgRNA info', 'strand')] == '-' else 0)
            if len(geneSgIds) < sgRNAsToPick and row['off-target filters'].loc[offTargetLevels[stringency]].all() \
                and matchREsites(oligoSeq, restrictionSites) \
                and checkOverlaps(leftPos, geneLeftPositions, nonoverlapMin):
                geneSgIds.append((sgId_v2,
                                  gene,transcript,
                                  row[('library table v2','sequence')], oligoSeq,
                                  row[('predicted score','predicted score')], np.nan,
                                 stringency))
                geneLeftPositions.append(leftPos)
                
        stringency += 1
            
    if len(geneSgIds) < sgRNAsToPick:
        unfinishedTss.append((gene, transcript)) #if the number of accepted sgRNAs is still less than sgRNAsToPick, discard gene
 #If there are not enough sgRNAs, report them all
    newSgIds.extend(geneSgIds)
