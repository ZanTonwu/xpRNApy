# Module xpRNApy
# Gene expression analysis of RNASeq data sets using genome annotation GTF and HTSeq-count read count data
# cool names include: Mustardo, Heizmappe, Xpression

import numpy as np
import HTSeq
import scipy.stats
# import scipy.interpolate.spline
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import pandas as pd
import seaborn
from matplotlib_venn import venn3, venn3_circles


class XpressionAnalysis(object):
    """
    This class represents an expression analysis project. Ggene annotation data in a dictionary to which expression data can be added.
    Several methods work on these data and include normalization functions as well as import and export tools. For convenience conversion
    functions to and from pandas DataFrames and CSV files are implemented. Heatmaps can be generated for various aspects of the dataset.

    Class variables:

        gene_annotation : list of str containing names of Keys and columns for gene annotation ["chrom", "start", "end", "strand", "length", "exon_length"]
        RBGcmap         : a matplotlib color map for heat maps red - black -green
        
    Data is held in the instance variables:

        genes (dict)     :  Uses gene_id to associate one dictionary per gene that contains annotation "chrom", "start", "end", "strand",
                            "length", "exon_length", as well as expression data associated with filenames
        filenames (list) :  Contains filenames of expression data files (HTSeq-count read count filed) added to the genes directorys

    Settings for controling behavious by instance variables:
    
        printVerboseLog (boolean) : Controls error and warning messages printed to the screen (default = True)
        
    """
    RBGcmap=LinearSegmentedColormap.from_list("mycmap",['#ff0000','#000000','#00ff00'])
    gene_annotation=["chrom", "start", "end", "strand", "length", "exon_length"]

    def __init__(self, GTF_filename=None, filenames=None, end_included=True): # end included refers to the genomic position of exon intervalls if number of iv.end includes the last base of the exon
        """
        This is the basic constructor for an expression analysis. It contains gene annotation data in a dictionary to which expression data can be added.
        Several methods work on these data and include normalization functions as well as import and export tools. For convenience conversion
        functions to and from pandas DataFrames and CSV files are implemented. Heatmaps can be generated for various aspects of the dataset.

        Parameters:

            GTF_filename (str)     : A valid filename for a GTF file that contains annotation for the genome of the analysis. This will be used to
                                     "chrom", "start", "end", "strand", "length", "exon_length" in the genes dictionary
            filenames (list)       : Contains filenames of expression data files (HTSeq-count read count filed) to be added to the genes directorys
            end_included (boolean) : Used by HTSeq to adjust the genomic position of exon intervalls from the GTF file if number of iv.end includes the last base of the exon
            
        """
        self.genes = dict()                              # dictionary indexed with gene_id containing dictionaries with gene data
        self.filenames = list()                          # list of filenames (str) containing read count data, used to identify datasets
        self.printVerboseLog=True                        # used to control screen print for error and warnings
        if not (GTF_filename is None):
            self.readGTFannotation(GTF_filename, end_included)
            if not (filenames is None):
                self.addReadCountDataFiles(filenames)

    def readGTFannotation(self, GTF_filename, end_included=True):
        """
        Read gene annotation data from a GTF file into a dictionary to which expression data can be added subsequently.

        Parameters:

            GTF_filename (str)     : A valid filename for a GTF file that contains annotation for the genome of the analysis. This will be used to
                                     "chrom", "start", "end", "strand", "length", "exon_length" in the genes dictionary
            end_included (boolean) : Used by HTSeq to adjust the genomic position of exon intervalls from the GTF file if number of iv.end includes the last base of the exon
            
        """
        GTF_file = HTSeq.GFF_Reader(GTF_filename, end_included)
        exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)     # used to store exon information as gene_id over genomic position
        for feature in GTF_file:
            if feature.type=='exon':
                gene_id = feature.attr['gene_id']
                chrom = feature.iv.chrom
                start = feature.iv.start
                end = feature.iv.end
                if start>end:
                    start, end = end, start
                strand = feature.iv.strand
                if gene_id in self.genes:
                    if self.genes[gene_id]['chrom']!=chrom:
                        self.genes[gene_id]['chrom']='ambiguous'
                        if self.printVerboseLog:
                            print "WARNING", gene_id,"has ambiguous chromosome"
                    else:
                        if self.genes[gene_id]['strand']==strand:
                            if self.genes[gene_id]['start']>start:
                                self.genes[gene_id]['start']=start
                                self.genes[gene_id]['length']=self.genes[gene_id]['end']-start
                            if self.genes[gene_id]['end']<end:
                                self.genes[gene_id]['end']=end
                                self.genes[gene_id]['length']=end-self.genes[gene_id]['start']
                            exons[feature.iv]+=gene_id
                        else:
                            if self.printVerboseLog:
                                print "WARNING",gene_id,"has ambiguous strand"
                else:
                    length=end-start
                    self.genes[gene_id]= {'chrom':chrom,'start':start, 'end':end, 'strand':strand, 'length':length}
                    exons[feature.iv]+=gene_id
        for iv, gene_set in exons.steps():                          # calculate total exon length per gene using the genomic array of sets exons
            length=iv.end-iv.start
            for gene_id in gene_set:
                if 'exon_length' in self.genes[gene_id]:
                    self.genes[gene_id]['exon_length']+=length
                else:
                    self.genes[gene_id]['exon_length']=length

    def removeReadCountData(self, filenames):
        """
        Removes read count data from the sample associated with filenames and removes the sample form filenames, but keeps all other read counts and
        gene annotation information.

        Parameters:

        filename (list) : Contains a list of filenames for which readcount data is to be removed
        """
        if type(filenames)==list:
            if len(filenames)>0:
                gene_info=XpressionAnalysis.gene_annotation
                for filename in filenames:
                    if filename not in gene_info:                     # check to not delete ["chrom", "start", "end", "strand", "length", "exon_length"]
                        for gene_id in self.genes:
                            self.genes[gene_id].pop(filename, None)   # delete the gene dict key if exist else silently pass None
        else:
            if self.printVerboseLog:
                print "WARNING removeReadCountData requires a list argument for filenames,", type(filenames), "was provided"

    def clearReadCountData(self):
        """
        Clears all read count data and sample filenames, but keeps gene annotation information.
        """
        newGenes = dict()
        gene_info=XpressionAnalysis.gene_annotation    # = ["chrom", "start", "end", "strand", "length", "exon_length"]
        for gene_id in self.genes:                     # make a copy of annotation information of current gene
            newGenes[gene_id]=dict()
            for info in gene_info:
                newGenes[gene_id][info]=self.genes[gene_id][info]
        self.genes=newGenes
        self.filenames=list()

    def clearAll(self):
        """
        Clears read count and gene annotation data as well as sample filenames from the analysis project. Similar to instanciating a new object.
        """
        self.genes = dict()
        self.filenames=list()

    def copy(self, filenames=None):
        """
        Copies gene annotation and read count data from the sample associated with filenames as well as filenames into a new XpressionAnalysis object.
        If no filenames are specified then all of the data are copied.

        Parameters:

        filename (list, None) : Contains a list of filenames for which readcount data is to be copied. If None (default) all of the data is copied.

        Returns:

            NewAnalysis (object) : New XpressionAnalysis with normalized read counts
            
        """
        newAnalysis=XpressionAnalysis()
        if filenames is None:
            for gene_id in self.genes:            # make a copy of the current gene
                newAnalysis.genes[gene_id]=dict()
                for info in self.genes[gene_id]:
                    newAnalysis.genes[gene_id][info]=self.genes[gene_id][info]            
            newAnalysis.filenames=self.filenames
        else:
            gene_info=XpressionAnalysis.gene_annotation # = ["chrom", "start", "end", "strand", "length", "exon_length"]
            gene_info=gene_info + filenames
            for gene_id in self.genes:            # make a copy of the annotation and count data for filenames of the current gene
                newAnalysis.genes[gene_id]=dict()
                for info in gene_info:
                    newAnalysis.genes[gene_id][info]=self.genes[gene_id][info]            
            newAnalysis.filenames=filenames
        return newAnalysis            

    def to_RPKM(self):
        """
        Converts read counts to reads per kilobase transcript using gene annotation exon_length, and per million reads in the dataset of the experiment (RPKM). This
        can be useful for normalizing datasets, but also see to_TPM() that might achieve a superior outcome.
        """
        totalA, totalX = self.getTotalReadCountsAX(self.filenames)  # get the readcounts for X and autosomes
        for gene_id in self.genes:
            if 'exon_length' in self.genes[gene_id]:
                samples=[i for i in self.genes[gene_id].keys() if i not in XpressionAnalysis.gene_annotation]  # get all data but not annotation keys from the gene dict
                for sample in samples:
                    rpk = float(self.genes[gene_id][sample]*1000) / float(self.genes[gene_id]['exon_length'])
                    rpkm = rpk * 1000000.0 / float(totalX[sample] + totalA[sample])
                    self.genes[gene_id][sample] = rpkm
            else:
                if self.printVerboseLog:
                    print "WARNING exon_length is not defined for",gene_id, ", RPKM could not be calculated"                

    def to_TPM(self):
        """
        Converts read counts to transcripts per million (TPM) using gene annotation exon_length. Transcripts per million are preferred over RPKM as gene read counts for
        genes are first normalized per kilobase of transcript. In a second step the data for one experiment is normalized to one million transcripts by dividing by the
        sum of reads per kilobase transcript counts for all genes and multiplying by 1000000. For normalization of expression data this can provide a more even weight
        of genes of vastly differing length - else large genes with high read counts likely dominate the normalization.
        """
        for gene_id in self.genes:
            if 'exon_length' in self.genes[gene_id]:
                samples=[i for i in self.genes[gene_id].keys() if i not in XpressionAnalysis.gene_annotation]  # get all data but not annotation keys from the gene dict
                for sample in samples:
                    rpk = float(self.genes[gene_id][sample]*1000) / float(self.genes[gene_id]['exon_length'])
                    self.genes[gene_id][sample] = rpk
                totalA, totalX = self.getTotalReadCountsAX(self.filenames)  # get the readcounts, which are already adjusted per kilobase of transcript, for X and autosomes
                for sample in samples:
                    rpk = float(self.genes[gene_id][sample])
                    tpm = rpk * 1000000.0 / float(totalX[sample] + totalA[sample])
                    self.genes[gene_id][sample] = tpm                    
            else:
                if self.printVerboseLog:
                    print "WARNING exon_length is not defined for",gene_id, ", TPM could not be calculated"                

    def addReadCountDataFiles(self, filenames, path="", suffix="_count"):
        """
        Add read count data from an HTSeq-count formated file to the current analysis. A path can be specified for a folder containing the count files, which is useful
        for keeping filenames succinct as these are used as column lables for identifying samples. A suffix is added to the filenale (default "_count") again to keep
        filenames meaninful sample identifiers.
        If a gene_id is present in the analysis the read count will be added using the filename as the key, if the key does not exist. If the key already exists or
        if the gene_id does not exist the read counts are not discarded - a warning is printed.

        Parameters:

            filenames (list) : Contains a list of filenames (str) with the read count data. The files contain lines formated <gene_id> <count> from HTSeq-count.
            
        """
        for filename in filenames:
            data_added=False
            if self.printVerboseLog:
                print "INFO Adding data from", filename
            with open(path + filename + suffix,"r") as count_file:
                for line in count_file:
                    gene_id,count=line.split()                      # each line contains exactly two items: gene_id and integer count data
                    if gene_id in self.genes:
                        if filename in self.genes[gene_id]:
                            if self.printVerboseLog:
                                print "WARNING",gene_id,"has existing",filename
                        else:
                            self.genes[gene_id][filename]= int(count)
                            data_added=True
                    else:
                        if self.printVerboseLog:
                            print "WARNING",gene_id, "not in genome"
            if data_added==True:
                self.filenames.append(filename)

    def addReadCountData(self, XprAnalysis, filenames=None, replaceData=False, mergeGenes=False):
        """
        Add data from an existing analysis object. The gene_ids are used for identification. If replace is True and filenames of the samples exists it will be
        replaced with new data, if replase is False the data remains unchanged, a warning is generated. This is useful for combining analyses of samples that
        have been either obtained from different experiments or processed separately. If analysis with different genes are to be merged mergeGenes can be set
        to True to combine genes, but care must be taken as different count data might then exist for different genes in the analysis, ie some filenames might
        not be defined for all genes of the resulting analysis dataset.

        Parameters:

        XprAnalysis (object)  : A XpressionAnalysis object from which data will be added to the current analysis
        filename (list, None) : Contains a list of filenames for which readcount data is to be added. If None (default) all of the data is copied.
        replace (boolean)     : Set to True to overwrite existing data for a gene with the new data to be added. Note that if filenames contains identifiers
                                for gene annotation this will also be overwritten. If False (default) existing data will not be overwritten, the new data will
                                be ignored, and a warning printed.
        mergeGenes (boolean)  : If False only data for genes that are present in the current analysis will be added. If True then genes that do not exist in
                                the current analysis will be added from the new data. This can be useful to add additional genes, but care must be taken as in
                                the resulting analysis dataset filenames might not be defined for all genes.
        
        """
        if filenames is None:
            filenameList = XprAnalysis.filenames
        else:
            filenameList = filenames
        if len(filenameList)>0:
            dataAdded = False
            for gene_id in XprAnalysis.genes:            # make a copy of the current gene
                if gene_id in self.genes:
                    for info in filenameList:
                        if (info not in self.genes[gene_id]) or replaceData:
                            self.genes[gene_id][info]=XprAnalysis.genes[gene_id][info]
                            dataAdded=True
                        else:
                            if printVerboseLog:
                                print "WARNING gene_id",gene_id,"has existing", info,"count data. New data will be skipped (replace=False)."
                else:
                    if mergeGenes==True:
                        self.genes[gene_id] = dict()    # add a new gene
                        gene_info=XpressionAnalysis.gene_annotation # = ["chrom", "start", "end", "strand", "length", "exon_length"]
                        gene_info=gene_info + filenameList
                        for info in gene_info:
                            self.genes[gene_id][info]=XprAnalysis.genes[gene_id][info]
                        dataAdded=True
                    else:
                        if printVerboseLog:
                            print "WARNING gene_id",gene_id,"is not present in gene annotation of current analysis. Gene will be skipped."
            if dataAdded==True:
                newFilenames=[fnm for fnm in fileNameList if fnm not in self.filenames]
                self.filenames = self.filenames + newFilenames
        else:
            if printVerboseLog:
                print "WARNING no read count data was found."
            

    def writeDataFile(self, filename, sep='\t', overwrite_existing_file=False):
        """
        Writes records from the genes directory including annotation info and read count data for datasets specified in filename to a csv file.

        Parameters:

            filename (str)                    : Filename of the output file. Must be a valid path. 
            sep (str)                         : Separator used for data fields (default = '\t')
            overwrite_existing_file (boolean) : If True an existing file will be overwritten (default = False)
            
        """
        if overwrite_existing_file==False:   # check if filename already exists and exit before overwriting
            try:
                with open(filename, "r") as inFile:
                    return
            except IOError:
                pass
        with open(filename,"w") as outfile:
            if self.printVerboseLog:
                print "Writing data to disk ..."
            columns=["gene_id"]
            columns=columns+XpressionAnalysis.gene_annotation # =["chrom", "start", "end", "strand", "length", "exon_length"]
            columns=columns+self.filenames
            header=sep.join(columns)
            header=header+"\n"
            outfile.write(header)
            num=0
            for gene_id in sorted(self.genes.keys()):
                row=[gene_id]
                for column in columns[1:]:
                    if column in self.genes[gene_id]:
                        row.append(str(self.genes[gene_id][column]))
                    else:
                        if self.printVerboseLog:
                            print "WARNING", filename, "is not defined for", gene_id
                        return
                record=sep.join(row)
                record=record+"\n"
                outfile.write(record)
                num=num+1
            if self.printVerboseLog:
                print num,"of records written to", filename

    def print_count_data(self, limit=20, chrom=None):
        """
        Prints a number of gene records from the genes directory including chrom annotation and read count data for datasets specified in filename.

        Parameters:

            limit (int, None) : Limit of number of genes to print records for. If limit is 0 or None all records will be printed
            chrom (str, None) : Print records for specified chromosome. Default is None prints for all chromosomes.
            
        """
        i=1
        for gene_id in sorted(self.genes.keys()):
            if not (chrom is None):
                if chrom!=self.genes[gene_id]['chrom']:
                    continue
            print i, gene_id, self.genes[gene_id]['chrom'],
            for filename in self.filenames:
                print filename,":", self.genes[gene_id][filename],
            print
            i=i+1
            if not (limit is None):
                if limit>0 and i>limit:
                    break

    def getTotalReadCountsAX(self, filenames):
        """
        Calculates the total number of reads in each dataset referenced by filenames and returns two dictionaries with counts of autosomal and X-chromosomal (chrX) read for each dataset.
        Keys are filenames containing the data sets.

        Parameters:

            filenames (list) : List of filenames (str) that reference datasets. All must be present in the genes directory. If any of the annotation Keys are used behaviour is undefined.

        Returns:

            totalA (dict) : Total autosmal read counts with Keys from filenames containing datasets
            totalX (dict) : Total X-chromosomal read counts as identified with 'chrX'
            
        """
        totalA = dict()
        totalX = dict()
        for filename in filenames:
            totalA[filename] = 0
            totalX[filename] = 0
        for gene_id in self.genes:
            for filename in filenames:
                if filename in self.genes[gene_id]:
                    if self.genes[gene_id]['chrom']=="chrX":
                        totalX[filename]+=self.genes[gene_id][filename]
                    else:
                        totalA[filename]+=self.genes[gene_id][filename]
                else:
                    if self.printVerboseLog:
                        print "WARNING", filename, "is not defined for", gene_id
                    return
        return totalA, totalX
    
    def getTotalChromReadCounts(self, filenames, chrom):
        """
        Calculates the total number of reads in each dataset referenced by filenames and returns a dictionary with counts of read mapping to the specified chromosome for each dataset.
        Keys are filenames containing the data sets.

        Parameters:

            filenames (list) : List of filenames (str) that reference datasets. All must be present in the genes directory. If any of the annotation Keys are used behaviour is undefined.

        Returns:

            total (dict) : Total chromosomal read counts with Keys from filenames containing datasets
            
        """
        total = dict()
        for filename in filenames:
            total[filename] = 0
        for gene_id in self.genes:
            for filename in filenames:
                if filename in self.genes[gene_id]:
                    if self.genes[gene_id]['chrom']==chrom:
                        total[filename]+=self.genes[gene_id][filename]
                else:
                    if self.printVerboseLog:
                        print "WARNING", filename, "is not defined for", gene_id
                    return
        return total

    def getLinRegSI(self, filenames, autosomal=True, zero_intercept=True, reference=0):
        """
        Calculate a linear regression for normalizing datasets in filenames using one dataset as reference. If zero_intercept is set to False then 
        slope, intercept, corr_coeff = scipy.stats.linregress(x,y=) is used. For normalizing the dataset contained in filename use the formula:

        # normalizedReadCount[gene_id] = ( readCount[gene_id] - intercept ) / slope

        Parameters:

            filenames (list)    : List of filenames (str) that reference datasets. All must be present in the genes directory. If any of the annotation Keys are used behaviour is undefined.
            autosomal (boolean) : Omit X-linked genes identified by chrom=='chrX' from linear regression (default = True)
            zero_intercept (boolean) : Calculate a linear regression line that goes through the origin (default = True). This is recommended.
            reference (int)     : Index to the dataset in filenames that should be used as a reference for calculating the linear regression

        Returns:

            slopes (dict)       : Slopes of the linear regression with Keys from filenames containing datasets
            intercepts (dict)   : Intercepts of the linear regression line

        """
        slopes = dict()                                         
        intercepts = dict()
        gene_expression = dict()     
        for filename in filenames:
            slopes[filename] = 0
            intercepts[filename] = 0
            gene_expression[filename] = list() # will hold read counts for genes for linear regression
        for gene_id in self.genes:
            if autosomal==False or self.genes[gene_id]['chrom']!="chrX": # omit x-linked genes if autosomal is true
                for filename in filenames:
                    if filename in self.genes[gene_id]:
                        gene_expression[filename].append(self.genes[gene_id][filename])
                    else:
                        if self.printVerboseLog:
                            print "WARNING", filename, "is not defined for", gene_id
                        return
        filename_reference = filenames[reference]
        for filename in filenames:
            if filename != filename_reference:
                if zero_intercept == True:
                    slope = self.lin_reg_slope_no_intercept(gene_expression[filename_reference], gene_expression[filename])
                    intercept = 0
                else:
                    slope, intercept, rval, pval, stderr = scipy.stats.linregress(gene_expression[filename_reference], y=gene_expression[filename])                
                slopes[filename] = slope
                intercepts[filename] = intercept
        return slopes, intercepts

    def getLinRegChrom(self, filenames, chrom, zero_intercept=True, reference=0):
        """
        Calculate a linear regression for genes on the specified chromosome for normalizing datasets in filenames using one dataset as reference.
        If zero_intercept is set to False then slope, intercept, corr_coeff = scipy.stats.linregress(x,y=) is used.
        For normalizing the dataset contained in filename use the formula:

        # normalizedReadCount[gene_id] = ( readCount[gene_id] - intercept ) / slope

        Parameters:

            filenames (list)         : List of filenames (str) that reference datasets. All must be present in the genes directory. If any of the
                                       annotation Keys are used behaviour is undefined.
            chrom (str)              : Only use genes from the specified chromosome eg. chrom=='chrX' for linear regression
            zero_intercept (boolean) : Calculate a linear regression line that goes through the origin (default = True). This is recommended.
            reference (int)          : Index to the dataset in filenames that should be used as a reference for calculating the linear regression

        Returns:

            slopes (dict)       : Slopes of the linear regression with Keys from filenames containing datasets
            intercepts (dict)   : Intercepts of the linear regression line

        """
        slopes = dict()                                         
        intercepts = dict()
        gene_expression = dict()     
        for filename in filenames:
            slopes[filename] = 0
            intercepts[filename] = 0
            gene_expression[filename] = list() # will hold read counts for genes for linear regression
        for gene_id in self.genes:
            if self.genes[gene_id]['chrom']==chrom: # use only genes on the specified chromosome for calculation
                for filename in filenames:
                    if filename in self.genes[gene_id]:
                        gene_expression[filename].append(self.genes[gene_id][filename])
                    else:
                        if self.printVerboseLog:
                            print "WARNING", filename, "is not defined for", gene_id
                        return
        filename_reference = filenames[reference]
        for filename in filenames:
            if filename != filename_reference:
                if zero_intercept == True:
                    slope = self.lin_reg_slope_no_intercept(gene_expression[filename_reference], gene_expression[filename])
                    intercept = 0
                else:
                    slope, intercept, rval, pval, stderr = scipy.stats.linregress(gene_expression[filename_reference], y=gene_expression[filename])                
                slopes[filename] = slope
                intercepts[filename] = intercept
        return slopes, intercepts

    def normalizedExpression(self, filenames, autosomal=True, LinReg=False, zero_intercept=True, reference=0, floorValue=0, floorCounts=True):
        """
        Returns a new XpressionAnalysis object with normalized read counts for all datasets specified in filenames. The gene annotation is also copied.
        Three normalization methods are supported: (1) Total read count, (2) Linear regression with zero intercept, (3) linear regression with intercept
        Normalalization parameters can be calculated from the entire list of genes or by omitting X-linked genes identified ('chrX'). For linear regression
        an index to the reference dataset can be supplied. Additionally, the expression values can be floored to floorValue. 0 or 1 could be useful to
        remove negative read counts and when attempting division or log transform.
        
        Parameters:

            filenames (list)         : List of filenames (str) that reference datasets. All must be present in the genes directory. If any of the annotation Keys
                                       are used behaviour is undefined.
            autosomal (boolean)      : Omit X-linked genes identified by chrom=='chrX' from calculation of normalization factors (default = True)
            LinReg (boolean)         : If False datasets are normalized by total read count method. Can be restricted to autosomal.
            zero_intercept (boolean) : Calculate a linear regression line that goes through the origin (default = True). This is recommended.
            reference (int)          : Index to the dataset in filenames that should be used as a reference for calculating the linear regression
            floorValue (int, float)  : Value to floor normalized read counts below floorValue (default = 0)
            floorCounts (boolean)    : If True normalized read counts will be floored to floorValue

        Returns:

            XprAnalysis (object) : New XpressionAnalysis with normalized read counts
        
        """
        normalizedGenes = dict()
        gene_info=XpressionAnalysis.gene_annotation # = ["chrom", "start", "end", "strand", "length", "exon_length"]
        for gene_id in self.genes:            # make a copy of annotation information of current gene
            normalizedGenes[gene_id]=dict()
            for info in gene_info:
                normalizedGenes[gene_id][info]=self.genes[gene_id][info]
        if LinReg==True:                      # calculate normalized read counts using linear regression methods
            slope, intercept = self.getLinRegSI(filenames, autosomal=autosomal, zero_intercept=zero_intercept, reference=reference)
            for gene_id in self.genes:
                for filename in filenames:
                    if intercept[filename]==0 and slope[filename]==0:
                        normalizedReadCount = float(self.genes[gene_id][filename]) # if slope and intercept 0,0 then REFERENCE copy values
                    else:
                        normalizedReadCount = float(self.genes[gene_id][filename] - intercept[filename]) / float(slope[filename])
                    if floorCounts==True and normalizedReadCount < float(floorValue):
                        normalizedReadCount = float(floorValue)         # avoid negative read counts by flooring to floor value, ie 0 or 1 can be useful
                    normalizedGenes[gene_id][filename] = normalizedReadCount
        else:                                 # calculate normalized read counts using total read count methods
            totalA, totalX = self.getTotalReadCountsAX(filenames)
            if autosomal==True:
                maximalReadCount = max([totalA[i] for i in totalA])
                for gene_id in self.genes:
                    for filename in filenames:
                        normalizedReadCount = float(self.genes[gene_id][filename] * maximalReadCount) / float(totalA[filename])
                        if floorCounts==True and normalizedReadCount < float(floorValue):
                            normalizedReadCount = float(floorValue)         # avoid negative read counts by flooring to floor value, ie 0 or 1 can be useful
                        normalizedGenes[gene_id][filename] = normalizedReadCount
            else:   # autosomal == False
                maximalReadCount = max([(totalA[i]+totalX[i]) for i in totalA])
                for gene_id in self.genes:
                    for filename in filenames:
                        normalizedReadCount = float(self.genes[gene_id][filename] * maximalReadCount) / float(totalA[filename]+totalX[filename])
                        if floorCounts==True and normalizedReadCount < float(floorValue):
                            normalizedReadCount = float(floorValue)         # avoid negative read counts by flooring to floor value, ie 0 or 1 can be useful
                        normalizedGenes[gene_id][filename] = normalizedReadCount
        normData = XpressionAnalysis()
        normData.genes=normalizedGenes
        normData.filenames=filenames
        return normData

    def normalizeChromosome(self, filenames, chrom, LinReg=False, zero_intercept=True, reference=0, ToCtrl=False, floorValue=0, floorCounts=True):
        """
        Returns a new XpressionAnalysis object with normalized read counts for genes on a specified chromosome for all datasets specified in filenames. The gene annotation is also copied.
        Three normalization methods are supported: (1) Total read count, (2) Linear regression with zero intercept, (3) linear regression with intercept
        Normalalization parameters can be calculated from the entire list of genes or by omitting X-linked genes identified ('chrX'). For linear regression
        an index to the reference dataset can be supplied. Additionally, the expression values can be floored to floorValue. 0 or 1 could be useful to
        remove negative read counts and when attempting division or log transform.
        Note: Normalization parameters are calculated from and applied to the genes of the specified chromosome, genes on other chromosomes are not changed. It is possible to provide an
              index to a reference samples through setting IdxToCtrl. This is useful if a cell line effect (chromosomal aneuploidy) is to be normalized by using normalization factors
              from control conditions. In this case normalization factors are calculated for samples i = 0..(n/2-1) and applied to sample i and i+(n/2). Samples are ordered in filenames.
        
        Parameters:

            filenames (list)         : List of filenames (str) that reference datasets. All must be present in the genes directory. If any of the annotation Keys
                                       are used behaviour is undefined.
            chrom (str)              : Chromosome for which expression of genes will be normalized. Only for the genes will be normalization factors calculated and applied.
            LinReg (boolean)         : If False datasets are normalized by total read count method. Can be restricted to autosomal.
            zero_intercept (boolean) : Calculate a linear regression line that goes through the origin (default = True). This is recommended.
            reference (int)          : Index to the dataset in filenames that should be used as a reference for calculating the linear regression
            ToCtrl (boolean)         : if provided the normalization factor will calculated for sample i will also be applied to sample i+(n/2)
                                       where n is the number of samples (filenames). Useful for normalizing sample cell line effects with value from control condition.
            floorValue (int, float)  : Value to floor normalized read counts below floorValue (default = 0)
            floorCounts (boolean)    : If True normalized read counts will be floored to floorValue

        Returns:

            XprAnalysis (object) : New XpressionAnalysis with normalized read counts for the specified chromosome
        
        """
        normalizedGenes = dict()
        gene_info=XpressionAnalysis.gene_annotation # = ["chrom", "start", "end", "strand", "length", "exon_length"]
        for gene_id in self.genes:            # make a copy of annotation information of current gene
            normalizedGenes[gene_id]=dict()
            for info in gene_info:
                normalizedGenes[gene_id][info]=self.genes[gene_id][info]
        if LinReg==True:                      # calculate normalized read counts using linear regression methods
            if ToCtrl == False:
                slope, intercept = self.getLinRegChrom(filenames, chrom, zero_intercept=zero_intercept, reference=reference)
            else:
                SampleNum = len(filenames)
                if SampleNum % 2 != 0:
                    if self.printVerboseLog:
                        print "normalizeChrom() WARNING: cannot use", SampleNum, "samples with normalization ToCtrl=True, requires even number of samples"
                    return None
                slope, intercept = self.getLinRegChrom(filenames[:SampleNum/2], chrom, zero_intercept=zero_intercept, reference=reference)
                for idx in range(SampleNum/2):
                    slope[filenames[idx+(SampleNum/2)]]=slope[filenames[idx]]
                    intercept[filenames[idx+(SampleNum/2)]]=slope[filenames[idx]]                
            for gene_id in self.genes:
                for filename in filenames:
                    if (intercept[filename]==0 and slope[filename]==0) or self.genes[gene_id]['chrom']!=chrom:
                        normalizedReadCount = float(self.genes[gene_id][filename]) # if slope and intercept 0,0 then REFERENCE or gene not on chromosome -> copy value
                    else:
                        normalizedReadCount = float(self.genes[gene_id][filename] - intercept[filename]) / float(slope[filename])
                    if floorCounts==True and normalizedReadCount < float(floorValue):
                        normalizedReadCount = float(floorValue)         # avoid negative read counts by flooring to floor value, ie 0 or 1 can be useful
                    normalizedGenes[gene_id][filename] = normalizedReadCount
        else:                                 # calculate normalized read counts using total read count methods
            if ToCtrl == False:
                total = self.getTotalChromReadCounts(filenames, chrom)
            else:
                SampleNum = len(filenames)
                if SampleNum % 2 != 0:
                    if self.printVerboseLog:
                        print "normalizeChrom() WARNING: cannot use", SampleNum, "samples with normalization ToCtrl=True, requires even number of samples"
                    return None
                total = self.getTotalChromReadCounts(filenames[:SampleNum/2], chrom)
                for idx in range(SampleNum/2):
                    total[filenames[idx+(SampleNum/2)]]=total[filenames[idx]]                
            maximalReadCount = max([total[i] for i in total])
            for gene_id in self.genes:
                for filename in filenames:
                    if self.genes[gene_id]['chrom']==chrom:
                        normalizedReadCount = float(self.genes[gene_id][filename] * maximalReadCount) / float(total[filename])
                    else:
                        normalizedReadCount = float(self.genes[gene_id][filename]) # if gene not on chromosome -> copy value
                    if floorCounts==True and normalizedReadCount < float(floorValue):
                        normalizedReadCount = float(floorValue)         # avoid negative read counts by flooring to floor value, ie 0 or 1 can be useful
                    normalizedGenes[gene_id][filename] = normalizedReadCount
        normData = XpressionAnalysis()
        normData.genes=normalizedGenes
        normData.filenames=filenames
        return normData

    def Xpression(self, filenames, chrom="chrX", sortGenePos=True):
        """
        Returns gene expression data for one chromosome in a numpy Array. Genes are in rows and can be sorted by chromosomal position.
        Data in columns correspond to datasets referenced by filenames. Expression values can be read counts or normalized readcounts of type float.
        Index as follows:

        # Array[gene_id, filename]

        Parameters:

            filenames (list)       : List of filenames (str) that reference datasets. All must be present in the genes dictionary.
                                     If any of the annotation Keys are used behaviour is undefined.
            chrom (str)            : Specified the chromosome for which expression data is generated (default chrom='chrX')
            sortGenePos (boolean)  : Sort the expression data according to the position of the genes on the chromosome
            reference (int)        : Index to the dataset in filenames that should be used as a reference for calculating the linear regression

        Returns:

            exprArray (numpy.Array) : A numpy Array with gene expression data. Genes are in rows, datasets correspond to columns.

        """
        geneList=[(i, self.genes[i]['start']) for i in self.genes if self.genes[i]['chrom']==chrom]
        if sortGenePos==True:
            geneList=sorted(geneList, key=lambda x: x[1])
        columns=len(filenames)
        rows=len(geneList)
        if (rows<=0) or (columns<=0):
            if self.printVerboseLog:
                print "Xpression()  WARNING no expression data found for",chrom
            return None
        exprArray = np.zeros([rows,columns],dtype=float)
        ridx=0
        for gene_id, pos in geneList:
            cidx=0
            for filename in filenames:
                if filename in self.genes[gene_id]:
                    exprArray[ridx,cidx]=self.genes[gene_id][filename]
                    cidx+=1
                else:
                    if self.printVerboseLog:
                        print "Xpression() WARNING",filename,"not defined for", gene_id
                    return None
            ridx+=1
        return exprArray

    def get_filenames_from_dict(self):
        for gene_id in self.genes:
            fnames=list(self.genes[gene_id].keys())
            break                           # we only need access one gene_id to get the keys corresponding to column headings
        for info in self.gene_annotation:   # now subtract the keys corresponding to annotation and not data
            fnames.remove(info)
        return fnames

    def readFilenames(self, filename):
        fnames=list()
        with open(filename, "r") as f:
            for l in f:
                if len(l.rstrip())>0:
                    fnames.append(l.rstrip())
        self.filenames=fnames

    def writeFilenames(self, filename):
        with open(filename, "w") as outfile:
            for fname in self.filenames:
                outfile.write(fname+"\n")

    def getDataFrame(self):
        columns=XpressionAnalysis.gene_annotation # =["chrom", "start", "end", "strand", "length", "exon_length"]
        columns=columns+self.filenames
        DF=pd.DataFrame().from_dict(self.genes, orient='index')
        DF=DF[columns]
        DF.index.name="gene_id"
        return DF

    def setGenesDictFromDataFrame(self, df):
        self.genes = df.to_dict(orient="index")
        fnames = list(df.columns.values.tolist())
        self.filenames = fnames
        for info in self.gene_annotation:   # now subtract the keys corresponding to annotation and not data
            fnames.remove(info)
        if "gene_id" in fnames:
            fnames.remove("gene_id")
        self.filenames=fnames

    def readGenesDictFromCSV(self, CSV_filename, sep='\t'):
        df = pd.DataFrame().from_csv(CSV_filename, sep=sep)
        self.setGenesDictFromDataFrame(df)

    @classmethod
    def from_dict(cls, genes_dict):
        XprAnalysis = cls()
        XprAnalysis.genes=genes_dict
        XprAnalysis.filenames = XprAnalysis.get_filenames_from_dict()
        return XprAnalysis

    @classmethod
    def from_csv(cls, CSV_filename, sep='\t'):
        XprAnalysis = cls()
        XprAnalysis.readGenesDictFromCSV(CSV_filename, sep)
        return XprAnalysis

    @classmethod
    def from_df(cls, df):
        XprAnalysis = cls()
        XprAnalysis.setGenesDictFromDataFrame(df)
        return XprAnalysis

    @staticmethod
    def lin_reg_slope_no_intercept(x_val, y_val):
        """
        Returns the slope of a regression line for series in x_val and y_val that goes through the origin [y = m*x]
        """
        x = np.array(x_val)
        y = np.array(y_val)
        return (np.mean(x*y)/np.mean(x*x))

    @staticmethod
    def normalizedXpression(exprArry, LinReg=True):
        """
        Uses sample pairs to adjust X-linked gene expression by estimating the number of Xs from no dox expression cols [0..(n-1)] relative to reference in col 0
        to correct the no dox and dox treated sample [n+sampleNum] with the no dox adjustment. Data should be normalized read counts.
        Requires sample pairs in column n (no dox) and column n + sampleNum (dox).

        Parameters:

            exprArry (numpy Array) : A numpy Array with gene expression data. Genes are in rows, datasets correspond to columns. Returned by Xpression()

        Returns:

            Array (numpy Array)     : With X-chromosomal normalized exression data in same format as input Array

        """
        (rows, cols) = exprArry.shape
        samples = cols/2
        if cols>samples*2:
            print "WARNING number of files",cols,"is not even. Requires",samples,"sample pairs for normalizing Xs"
            return
        slope=np.zeros([cols], dtype=float)
        for sample in range(samples):
            if LinReg==True:
                if sample>0:
                    m=XpressionAnalysis.lin_reg_slope_no_intercept(exprArry[:,0], exprArry[:,sample])
                else:
                    m=1
            else:
                m=np.sum(exprArry[:,sample])
            slope[sample]=m
            slope[sample+samples]=m
        if LinReg==True:
            return (exprArry / slope)
        else:
            return (exprArry * (max(slope)/slope))

    def HeatmapChroms(self, filenames, chromList, filename=None):
        numChrom=len(chromList)
        if numChrom<=0:
            if self.printVerboseLog:
                print "chromHeatmaps() WARNING nothing to plot: need chromosomes, <",chromList,"> specified"
            return
        fig, ax = plt.subplots(1, numChrom+1, gridspec_kw = {'width_ratios':[18]*numChrom + [1]}, sharex=False, sharey=False)
        for idx in range(numChrom):
            ax[idx].set_title(chromList[idx], loc='left', y=1.08)
            ax[idx].title.set_fontsize(20)
            exprA=self.Xpression(filenames, chrom=chromList[idx])
            zmapA=scipy.stats.zscore(exprA, axis=1)
            zmapA=zmapA[~np.isnan(zmapA).any(axis=1)]
            if idx<(numChrom-1):
                seaborn.heatmap(zmapA, cmap=self.RBGcmap, xticklabels=False, yticklabels=False, ax=ax[idx], cbar=False)
            else: # last set cbar
                seaborn.heatmap(zmapA, cmap=self.RBGcmap, xticklabels=False, yticklabels=False, ax=ax[idx], cbar=True, cbar_ax=ax[numChrom])
            ax[idx].xaxis.set_tick_params(labeltop='on', labelbottom='off')
            ax[idx].set_xticklabels(filenames)
            for label in ax[idx].get_xticklabels():
                label.set_rotation(90)
                label.set_fontsize(12)
        if (filename is None):
            plt.show()
        else:
            fig.savefig(filename)


    def HeizMappe(self, filenames, LinReg=True, chrom="chr1", filename=None):
        """
        Generated a pyplot figure using seaborn for X-chromosomal gene expression and a reference autosome. If a filename is given the figure is saved else displayed
        using matplotlib.

        Parameters:

            LinReg (boolean)     : Specified to use linear regression to adjust for X-linked gene expression between the samples using the no dox samples for estimating
                                   the number of X-chromosomes for both dox and no dox datasets of one sample pair
            chrom (str)          : The reference chromosome to display as a control (default = chr1)
            filename (str, None) : A valid file path to save the figure to. If None the figure will be displayed
            
        """
        exprX=self.Xpression(filenames)
        zmapX=scipy.stats.zscore(exprX, axis=1)
        zmapX=zmapX[~np.isnan(zmapX).any(axis=1)]
        exprA=self.Xpression(filenames, chrom=chrom)
        zmapA=scipy.stats.zscore(exprA, axis=1)
        zmapA=zmapA[~np.isnan(zmapA).any(axis=1)]
        exprXnorm=self.normalizedXpression(exprX, LinReg=LinReg)
        zmapXnorm=scipy.stats.zscore(exprXnorm, axis=1)
        zmapXnorm=zmapXnorm[~np.isnan(zmapXnorm).any(axis=1)]

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, gridspec_kw = {'width_ratios':[18,18,18,1]}, sharex=False, sharey=False)
        ax1.set_title("chrX", loc='left', y=1.08)
        ax1.title.set_fontsize(20)
        seaborn.heatmap(zmapX, cmap=self.RBGcmap, xticklabels=False, yticklabels=False, ax=ax1, cbar=False)
        ax1.xaxis.set_tick_params(labeltop='on', labelbottom='off')
        ax1.set_xticklabels(filenames)
        for label in ax1.get_xticklabels():
            label.set_rotation(90)
            label.set_fontsize(12)
        ax2.set_title("chr1", loc='left', y=1.08)
        ax2.title.set_fontsize(20)
        seaborn.heatmap(zmapA, cmap=self.RBGcmap, xticklabels=False, yticklabels=False, ax=ax2, cbar=False)
        ax2.xaxis.set_tick_params(labeltop='on', labelbottom='off')
        ax2.set_xticklabels(filenames)
        for label in ax2.get_xticklabels():
            label.set_rotation(90)
            label.set_fontsize(12)
        ax3.set_title("chrX normalized Xs", loc='left', y=1.08)
        ax3.title.set_fontsize(20)
        seaborn.heatmap(zmapXnorm, cmap=self.RBGcmap, xticklabels=False, yticklabels=False, ax=ax3, cbar=True, cbar_ax=ax4)
        ax3.xaxis.set_tick_params(labeltop='on', labelbottom='off')
        ax3.set_xticklabels(filenames)
        for label in ax3.get_xticklabels():
            label.set_rotation(90)
            label.set_fontsize(12)
        if (filename is None):
            plt.show()
        else:
            fig.savefig(filename)

# Analysis of RNAseq datasets associated with Monfort et al. --------------------------------------------------------------------------------------

def generateSampleGroups():
    samples={}
    for genotype in genotypes:
        samples[genotype] = dict()
        for condition in conditions:
            if condition == "Dox":
                dox_samples=list()
                for cellLine in cellLines[genotype]:
                    dox_samples.append(cellLine + "_DOX")
                samples[genotype][condition] = dox_samples
            else:
                samples[genotype][condition] = cellLines[genotype]
    return samples

def printSampleGroups():
    print "%10s   %40s   %40s" % ( "Genotype".ljust(10), conditions[0].ljust(40), conditions[1].ljust(40) )
    print "-"*(10+40+40)
    for genotype in genotypes:
        print "%10s   %40s   %40s" % ( genotype.ljust(10), " ; ".join(samples[genotype][conditions[0]]).ljust(40), " ; ".join(samples[genotype][conditions[1]]).ljust(40) )

def getSampleList(condition=None):
    sampleList=[]
    for genotype in genotypes:
        if condition is None:
            for c in conditions:
                sampleList=sampleList + samples[genotype][c]
        else:
            sampleList=sampleList + samples[genotype][condition]
    return sampleList

def getPairedSampleList():
    sampleList=[]
    for condition in conditions:
        for genotype in genotypes:            
            sampleList=sampleList + samples[genotype][condition]
    return sampleList

def getCountFilename(sample):
    return sample+"_counts"

def getCountFileList(folder=""):
    countFileList=[]
    for sample in getSampleList():
        countFileList.append(folder+getCountFilename(sample))
    return countFileList
    
def getTPM_DF(analysis):
    """
    Returns a pandas DataFrame with the expression data normalized to transcripts per million (TPM). The function accepts either an XpressionAnalysis object or a pandas DataFrame, whose columns
    must be initialized accordingly. Use getDataFrame to obtaina DataFrame from an XpressionAnalysis object. The annotation "exon_length" is used to noramlize to per kilobase transcript length.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    tpmDF=pd.DataFrame(columns=analysisDF.columns)
    for a in XpressionAnalysis.gene_annotation:         # copy annotation columns
        tpmDF[a]=analysisDF[a]
    for s in getSampleList():                           # calculate reads per kilobase of transcript (exon_length) = RPK
        tpmDF[s]=analysisDF[s]*1000/analysisDF["exon_length"]
    for s in getSampleList():                           # divide by sum of RPK for all genes of experiment time 1 mio to get TPM
        tpmDF[s]=tpmDF[s] * 1000000 / tpmDF[s].sum()
    return tpmDF

def getFold_DF(analysis, zero_floor_value=None):
    """
    Returns a pandas DataFrame with the fold difference for samples between the two conditions [0] and [1]. The function accepts either an XpressionAnalysis object or a pandas DataFrame, whose columns
    must be initialized accordingly. Use getDataFrame to obtaina DataFrame from an XpressionAnalysis object.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    cols=XpressionAnalysis.gene_annotation
    for genotype in genotypes:
        cols=cols+samples[genotype][conditions[0]]
    foldDF=pd.DataFrame(columns=cols)
    for a in XpressionAnalysis.gene_annotation:         # copy annotation columns
        foldDF[a]=analysisDF[a]
    for genotype in genotypes:                          # calculate the fold difference for samples between the two conditions [0] and [1]
        for (s0, s1) in zip(samples[genotype][conditions[0]],samples[genotype][conditions[1]]):
            if zero_floor_value is None:
                foldDF[s0]=analysisDF[s1]/analysisDF[s0]
            else:
                foldDF[s0]=analysisDF[s1]/np.maximum(zero_floor_value, analysisDF[s0])                
    return foldDF

def getFilterMinExpression_DF(analysis, condition=None, minValue=0):
    """
    Filter for rows that contain genes, whose expression is equal or above minValue. The requirement can be restricted to a condition by providing either
    "Ctrl" or "Dox" as argument. Default is None requiring minValue in all samples.  It accepts a analysis DataFrame or XpressionAnalysis object as input.
    Returns a pandas DataFrame with the remaining rows. Useful for filtering out genes with low readcount or TPM.
    """

    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    sampleList=getSampleList(condition=condition)
    return analysisDF[ analysisDF[sampleList].ge(minValue).all(axis=1) ]    # boolean selection of rows in which all sample columns are greater or equal to minValue

def getFilterNaN_DF(analysis):
    """
    Filter for rows that contain any NaN, inf, and -inf values, which cannot be processed further.
    Returns a pandas DataFrame with the remaining rows.  It accepts a analysis DataFrame or XpressionAnalysis object as input.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    analysisDF=analysisDF.replace([np.inf, -np.inf], np.nan)
    return analysisDF.dropna(axis=0, how="any")

def foldHistogram(analysis, bins=10, smooth=True, density=True, show_median_lines=True, filename="XA_fold_histogram.pdf"):
    """
    Generates a histogram of the fold values of genes for the chrX and autosomes. A figure with separate panels for each genotype is produced.
    Each panel contains density curves for the histograms of all replicates. It accepts a analysis DataFrame or XpressionAnalysis object as input.
    Median lines for fold Xchr and autosomal genes can be shown, the density curve can be smoothened by cubic polynomial interpolation. The number
    of bins for generating the histogram can be chosen. Note: Values outside -0.5 and +1.5 will be segregated as outliers in cutoff bins. X axis is
    trimmed to the range 0 to 1.2 as this is the interesting interval for repression by Xist (negative values are artefacts, and higher expression
    than 1 is expected to be a secondary effect but not related to Xist repression directly, for instance cell death related genes could be cases.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    fold_X_DF=analysisDF[analysisDF["chrom"]=="chrX"]
    fold_A_DF=analysisDF[analysisDF["chrom"]!="chrX"]
    fig=plt.figure(num="XA_fold_histogram", figsize=(4, 10), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=6, ncols=1, height_ratios=[1,1,1,1,1,1])
    colorList = ["darkred", "firebrick", "red", "steelblue", "seagreen", "darkcyan"]
    panel_row=0
    for genotype in genotypes:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 0])
            ax.set_title(genotype)
            c=0
            for s in samples[genotype]["Ctrl"]:
                y,x = np.histogram(fold_X_DF[s].values, density=density, bins=[-.6]+list(np.linspace(-.5,1.5,bins))+[2000])
                Xmedian_fold=fold_X_DF[s].median()
                Amedian_fold=fold_A_DF[s].median()
                x=x[1:-2]
                y=y[1:-1]
                if smooth:
                    x_smooth = np.linspace(x.min(), x.max(), 200)
                    y_smooth = scipy.interpolate.interp1d(x,y, kind='cubic') # scipy.interpolate.spline(x, y, x_smooth)
                    ax.fill_between(x_smooth, y_smooth(x_smooth), color=colorList[c], alpha=0.3)
                    ax.plot(x_smooth, y_smooth(x_smooth), color=colorList[c])
                else:
                    ax.fill_between(x, y, color=colorList[c], alpha=0.3)
                    ax.plot(x, y, color=colorList[c])
                y,x = np.histogram(fold_A_DF[s].values, density=density, bins=[-.6]+list(np.linspace(-.5,1.5,bins))+[2000])
                x=x[1:-2]
                y=y[1:-1]
                if smooth:
                    x_smooth = np.linspace(x.min(), x.max(), 200)
                    y_smooth = scipy.interpolate.interp1d(x,y, kind='cubic') # scipy.interpolate.spline(x, y, x_smooth)
                    ax.plot(x_smooth, y_smooth(x_smooth), color=colorList[c+3])
                else:
                    ax.plot(x, y, color=colorList[c+3])
                if show_median_lines:
                    ax.axvline(x=Xmedian_fold, ls="--", color=colorList[c])
                    ax.axvline(x=Amedian_fold, ls="--", color=colorList[c+3])
                c+=1
            ax.set_xlim(left=0, right=1.2)
            ax.set_ylim(bottom=0)
            ax.set_xticks([0,0.5,1])
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax.tick_params(which='minor', length=4, color='black')
            ax.xaxis.set_ticks_position('bottom')
            if show_median_lines==False:
                ax.axvline(x=1, alpha=0.5, ls="--", color="black")
        panel_row+=1
    fig.tight_layout()
    if smooth:
        fig.savefig(DATAPATH+"smooth_"+filename)
    else:
        fig.savefig(DATAPATH+filename)
    plt.close("XA_fold_histogram")

def averageFoldHistogram(analysis, bins=10, smooth=True, density=True, show_median_lines=True, filename="XA_averageFold_histogram.pdf"):
    """
    Generates a histogram of the fold values of genes for the chrX and autosomes. A figure with separate panels for each genotype is produced.
    Each panel contains density curves for the histograms of all replicates. It accepts a analysis DataFrame or XpressionAnalysis object as input.
    Median lines for fold Xchr and autosomal genes can be shown, the density curve can be smoothened by cubic polynomial interpolation. The number
    of bins for generating the histogram can be chosen. Note: Values outside -0.5 and +1.5 will be segregated as outliers in cutoff bins. X axis is
    trimmed to the range 0 to 1.2 as this is the interesting interval for repression by Xist (negative values are artefacts, and higher expression
    than 1 is expected to be a secondary effect but not related to Xist repression directly, for instance cell death related genes could be cases.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    fold_X_DF=analysisDF[analysisDF["chrom"]=="chrX"]
    fold_A_DF=analysisDF[analysisDF["chrom"]!="chrX"]
    fig=plt.figure(num="XA_fold_histogram", figsize=(4, 10), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=6, ncols=1, height_ratios=[1,1,1,1,1,1])
    colorList = ["darkred", "firebrick", "red", "steelblue", "seagreen", "darkcyan"]
    panel_row=0
    for genotype in genotypes:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 0])
            ax.set_title(genotype)
            c=0
            s_cols=[]
            for s in samples[genotype]["Ctrl"]:
                s_cols.append(s)
            mean_fold_X=fold_X_DF[s_cols].mean(axis=1)
            mean_fold_A=fold_A_DF[s_cols].mean(axis=1)
            y,x = np.histogram(mean_fold_X.values, density=density, bins=[-.6]+list(np.linspace(-.5,1.5,bins))+[2000])
            Xmedian_fold=mean_fold_X.median()
            Amedian_fold=mean_fold_A.median()
            x=x[1:-2]
            y=y[1:-1]
            if smooth:
                x_smooth = np.linspace(x.min(), x.max(), 200)
                y_smooth = scipy.interpolate.interp1d(x,y, kind='cubic') # scipy.interpolate.spline(x, y, x_smooth)
                #ax.fill_between(x_smooth, y_smooth(x_smooth), color=colorList[c], alpha=0.3)
                ax.plot(x_smooth, y_smooth(x_smooth), color=colorList[c])
            else:
                #ax.fill_between(x, y, color=colorList[c], alpha=0.3)
                ax.plot(x, y, color=colorList[c])
            y,x = np.histogram(mean_fold_A.values, density=density, bins=[-.6]+list(np.linspace(-.5,1.5,bins))+[2000])
            x=x[1:-2]
            y=y[1:-1]
            if smooth:
                x_smooth = np.linspace(x.min(), x.max(), 200)
                y_smooth = scipy.interpolate.interp1d(x,y, kind='cubic') # scipy.interpolate.spline(x, y, x_smooth)
                ax.plot(x_smooth, y_smooth(x_smooth), color=colorList[c+3])
            else:
                ax.plot(x, y, color=colorList[c+3])
            if show_median_lines:
                ax.axvline(x=Xmedian_fold, ls="--", color=colorList[c])
                ax.axvline(x=Amedian_fold, ls="--", color=colorList[c+3])
            ax.set_xlim(left=0, right=1.2)
            ax.set_ylim(bottom=0)
            ax.set_xticks([0,0.5,1])
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax.tick_params(which='minor', length=4, color='black')
            ax.xaxis.set_ticks_position('bottom')
            if show_median_lines==False:
                ax.axvline(x=1, alpha=0.5, ls="--", color="black")
        panel_row+=1
    fig.tight_layout()
    if smooth:
        fig.savefig(DATAPATH+"smooth_"+filename)
    else:
        fig.savefig(DATAPATH+filename)
    plt.close("XA_fold_histogram")

def getAverageExpressionLevel(analysis, calculate_median=False):
    """
    Calculates the mean or median expression level in all control samples (replicates).
    Returns a dataframe with index gene names and one column. It accepts a analysis DataFrame or XpressionAnalysis object as input.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    if calculate_median:
        return analysisDF[getSampleList(condition="Ctrl")].median(axis=1)
    else:
        return analysisDF[getSampleList(condition="Ctrl")].mean(axis=1)

def getAverageFold(analysis, calculate_median=False):
    """
    Calculates the mean or median fold for each genotype from replicates. It accepts a analysis DataFrame or XpressionAnalysis object as input.
    Returns a DataFrame with index gene names, columns "chrom" and "start", and each a column for each genotype with the mean or median fold.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    averageFoldDF=pd.DataFrame(columns=["chrom", "start"]+genotypes)
    averageFoldDF["chrom"]=analysisDF["chrom"]
    averageFoldDF["start"]=analysisDF["start"]
    foldDF=getFold_DF(analysisDF)
    for genotype in genotypes:
        if calculate_median:
            averageFoldDF[genotype]=foldDF[ samples[genotype][conditions[0]] ].median(axis=1)         # calculate the median fold for each genotype from replicates
        else:
            averageFoldDF[genotype]=foldDF[ samples[genotype][conditions[0]] ].mean(axis=1)           # calculate the mean fold for each genotype from replicates
    return averageFoldDF

def fexcorrPlot(analysis, show_median_lines=True, calculate_median=False, filename="fold_expression_correlation.pdf"):
    """
    Plots a correlation scatter plot of average fold (y) against average expression level per gene. It used logarithmic axes and lables Xist.
    It accepts a analysis DataFrame or XpressionAnalysis object as input.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    fexcorrDF=pd.DataFrame(columns=["chrom", "start", "expression_level"]+genotypes)
    fexcorrDF[["chrom","start"]+genotypes]=getAverageFold(analysisDF, calculate_median=calculate_median)
    fexcorrDF["expression_level"]=getAverageExpressionLevel(analysisDF, calculate_median=calculate_median) # calculate the average expression level in all control samples (replicates)
    fexcorrFilteredDF=getFilterNaN_DF(fexcorrDF)
    fexcorr_X_DF=fexcorrFilteredDF[fexcorrFilteredDF["chrom"]=="chrX"]
    fexcorr_A_DF=fexcorrFilteredDF[fexcorrFilteredDF["chrom"]!="chrX"]
    fig=plt.figure(num="fold_expression_correlation", figsize=(8, 12), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1,1,1], width_ratios=[1,1])
    panel_row=0
    panel_col=0
    for genotype in genotypes:
        if "DESPARATE_DEBUG" in globals():
            print genotype,"-----------------------" # try to help out in desparate situation
            print fexcorrFilteredDF.sort_values(by=[genotype], ascending=False).iloc[0:10]
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(genotype)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.plot(fexcorr_A_DF["expression_level"].values, fexcorr_A_DF[genotype].values, ".", color="darkcyan")
            ax.plot(fexcorr_X_DF["expression_level"].values, fexcorr_X_DF[genotype].values, ".", color="firebrick")
            ax.plot(fexcorr_X_DF.loc["Xist"]["expression_level"], fexcorr_X_DF.loc["Xist"][genotype], ".", color="black", zorder=3)
            ax.text(fexcorr_X_DF.loc["Xist"]["expression_level"]*1.5, fexcorr_X_DF.loc["Xist"][genotype]*0.85, "Xist",zorder=4)
            X_fold_median=fexcorr_X_DF[genotype].median(axis=0)
            A_fold_median=fexcorr_A_DF[genotype].median(axis=0)
            if show_median_lines==False:
                ax.axhline(y=1, alpha=0.5, ls="--", color="black")
            else:
                ax.axhline(y=X_fold_median, ls="--", color="firebrick")
                ax.axhline(y=A_fold_median, ls="--", color="darkcyan")
            ax.set_ylim(bottom=0.01, top=1000)
        panel_row+=1
        if panel_row>=3:
            panel_row=0
            panel_col+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("fold_expression_correlation")

def fexGroupBarPlot(analysis, calculate_median=False, filename="fold_expression_correlation.pdf"):
    """
    Plots a correlation scatter plot of average fold (y) against average expression level per gene. It used logarithmic axes and lables Xist.
    It accepts a analysis DataFrame or XpressionAnalysis object as input.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    fexcorrDF=pd.DataFrame(columns=["chrom", "start", "expression_level"]+genotypes)
    fexcorrDF[["chrom","start"]+genotypes]=getAverageFold(analysisDF, calculate_median=calculate_median)
    fexcorrDF["expression_level"]=getAverageExpressionLevel(analysisDF, calculate_median=calculate_median) # calculate the average expression level in all control samples (replicates)
    fexcorrFilteredDF=getFilterNaN_DF(fexcorrDF)
    grpLen=int(len(fexcorrFilteredDF.index)/3)
    fexcorrFilteredDF.sort_values("expression_level", axis = 0, ascending = False, inplace = True)
    high=fexcorrFilteredDF.iloc[0:grpLen]
    medium=fexcorrFilteredDF.iloc[grpLen:2*grpLen]
    low=fexcorrFilteredDF.iloc[2*grpLen:]
    fexcorr_X_DF=fexcorrFilteredDF[fexcorrFilteredDF["chrom"]=="chrX"]
    fexcorr_A_DF=fexcorrFilteredDF[fexcorrFilteredDF["chrom"]!="chrX"]
    fig=plt.figure(num="fold_expression_groups", figsize=(8, 12), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1,1,1], width_ratios=[1,1])
    panel_row=0
    panel_col=0
    for genotype in genotypes:
        Xmean=[ high[high["chrom"]=="chrX"][genotype].mean(), medium[medium["chrom"]=="chrX"][genotype].mean(), low[low["chrom"]=="chrX"][genotype].mean() ]
        Xstd=[ high[high["chrom"]=="chrX"][genotype].std(), medium[medium["chrom"]=="chrX"][genotype].std(), low[low["chrom"]=="chrX"][genotype].std() ]
        Amean=[ high[high["chrom"]!="chrX"][genotype].mean(), medium[medium["chrom"]!="chrX"][genotype].mean(), low[low["chrom"]!="chrX"][genotype].mean() ]
        Astd=[ high[high["chrom"]!="chrX"][genotype].std(), medium[medium["chrom"]!="chrX"][genotype].std(), low[low["chrom"]!="chrX"][genotype].std() ]
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(genotype)
            #ax.set_xscale("log")
            #ax.set_yscale("log")
            ax.bar([0,2,4], Xmean, yerr=Xstd, align='center', color="firebrick", ecolor='black', capsize=10)
            ax.bar([1,3,5], Amean, yerr=Astd, align='center', color="darkcyan", ecolor='black', capsize=10)
            ax.set_xticks([0.5,2.5,4.5])
            ax.set_xticklabels(["high", "medium", "low"])
            ax.set_ylabel('+Dox/Ctrl')
            ax.set_ylim(bottom=0, top=2.5)
        panel_row+=1
        if panel_row>=3:
            panel_row=0
            panel_col+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("fold_expression_groups")
    
def XposFoldPlot(analysis, show_median_lines=True, filename="fold_Xpos_plot.pdf"):
    """
    Plots the average fold (y) against X chromosomal position per gene. It used logarithmic y axes for fold, million unity for position (x), and lables Xist.
    It accepts a analysis DataFrame or XpressionAnalysis object as input.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    XposFoldDF=getAverageFold(analysisDF)
    XposFoldFilteredDF=getFilterNaN_DF(XposFoldDF)
    XposFold_X_DF=XposFoldFilteredDF[XposFoldFilteredDF["chrom"]=="chrX"]
    fig=plt.figure(num="fold_Xpos_plot", figsize=(8, 12), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1,1,1], width_ratios=[1,1])
    panel_row=0
    panel_col=0
    for genotype in genotypes:
        if "DESPARATE_DEBUG" in globals():
            print genotype,"-----------------------" # try to help out in desparate situation
            print XposFoldFilteredDF.sort_values(by=[genotype], ascending=False).iloc[0:10]
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(genotype)
            #ax.set_xscale("log")
            ax.set_yscale("log")
            plt.plot(XposFold_X_DF["start"].values, XposFold_X_DF[genotype].values, ".", color="firebrick")
            plt.plot(XposFold_X_DF.loc["Xist"]["start"], XposFold_X_DF.loc["Xist"][genotype], ".", color="black", zorder=3)
            plt.text(XposFold_X_DF.loc["Xist"]["start"]+5000000, XposFold_X_DF.loc["Xist"][genotype]*0.85, "Xist",zorder=4)
            X_fold_median=XposFold_X_DF[genotype].median(axis=0)
            ax.axhline(y=1, alpha=0.5, ls="--", color="black")
            if show_median_lines:
                ax.axhline(y=X_fold_median, ls="--", color="firebrick")
            ax.set_ylim(bottom=0.01, top=1000)
            xticks=list(range(0,max(XposFold_X_DF["start"].values),20000000))
            xtick_labels=list(range(0,max(xticks)/1000000+1,20))
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels)
            ax.tick_params(which='major', length=4, color='black')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        panel_row+=1
        if panel_row>=3:
            panel_row=0
            panel_col+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("fold_Xpos_plot")

def newXAnalysis():
    analysis = XpressionAnalysis(DATAPATH + "genes.gtf")
    analysis.addReadCountDataFiles(getSampleList(), path=DATAPATH, suffix="_counts")                 # use the sample identifiers as addReadCountDataFiles will construct the actual filenames
    analysis.writeDataFile(DATAPATH + "analysis_data", overwrite_existing_file=True)
    return analysis
    
if __name__ == "__main__":
    print "xpRNApy - X-linked gene eXPRession RNAseq Analysis in PYTHON"
    print "Usage:"
    print "\timport xpRNApy"
    print "\tanalysisProject = XpressionAnalysis.from_csv(\"data.txt\") ... to load an analysis"
    print "\tanalysisProject = newXAnalysis(DATAPATH=\"\\home\\linux\\Asun_FastQ\\\") ... to start a new analysis"

DATAPATH = "/home/linux/Asun_FastQ/" # "/home/linux/Asuns_Data/Asun_FastQ/" #
genotypes = ["WT", "Ubn1", "Ubn2", "Ubn1/2", "Hira", "Ubn1_2"]
conditions = ["Ctrl", "Dox"]
cellLines = { "WT": ["HATX3", "HATX8", "HATX12"],
              "Ubn1": ["U1_10", "U1_13", "U1_36"],
              "Ubn2": ["U9", "U14", "Ue38"],
              "Ubn1/2": ["U12_10", "U12_27", "U12_82"],
              "Ubn1_2": ["U12_8"],
              "Hira": ["H23", "HE46", "H18"] }
samples = generateSampleGroups()
analysisProject=newXAnalysis()
analysisDF=analysisProject.getDataFrame()
tpmDF=getTPM_DF(analysisDF)
foldDF=getFold_DF(tpmDF)
foldFilteredDF=getFilterNaN_DF(foldDF)

#foldHistogram(foldFilteredDF, smooth=False, show_median_lines=True, bins=20)
#foldHistogram(foldFilteredDF, smooth=True, show_median_lines=True, bins=20)

# analysis.HeizMappe(getPairedSampleList())
# tpm_analysis=analysis.copy()
# tpm_analysis.to_TPM()
# tpm_analysis.writeDataFile(DATAPATH + "analysis_TPM", overwrite_existing_file=True)
# tpm.analysis.HeizMappe(getPairedSampleList())

# fexcorrPlot(analysisProject)
#min10DF=getFilterMinExpression_DF(analysisProject, condition="Ctrl", minValue=10)
# fexcorrPlot(min10DF)
#min10tpmDF=getTPM_DF(min10DF)
#fexcorrPlot(min10tpmDF, show_median_lines=True)
#XposFoldPlot(tpmDF)


# Analysis of RNAseq and H3K27me3 ChIP correlations associated with Monfort et al. ------------------------------------------------------------------------

CHIP_DATA_PATH= "/home/linux/Asun_ChIPseq/Programs/" # "/home/linux/Asuns_Data/CHIPSeq_Asun2019/Programs/" #
ChIP_coverage_table="Analysis_f200r5q10_normAutCvg.csv"    # this is a csv file from the GeneCoverageAnalysis() function of a xChIPPy.XchipAnalysis object 
wt_ctrl_repl = ["WT_N8", "WT_N9", "WT_N10"]
wt_dox_repl = ["WT_DOX_N8", "WT_DOX_N9", "WT_DOX_N10"]
ko_ctrl_repl = ["KO_N8", "KO_N9", "KO_N10"]
ko_dox_repl = ["KO_DOX_N8", "KO_DOX_N9", "KO_DOX_N10"]
ChIP_samples = { "WT": { "Ctrl": wt_ctrl_repl, "Dox": wt_dox_repl },
                 "KO": { "Ctrl": ko_ctrl_repl, "Dox": ko_dox_repl } }
TSS_suffix = "AvgTSSCvg"
Gene_suffix = "AvgGeneCvg"
ChIP_DF=pd.DataFrame.from_csv(CHIP_DATA_PATH+ChIP_coverage_table)

def K27ExpFoldPlot(analysis, ChIP_DF, calculate_median=False, filename="K27_expr_fold.png"):
    """
    Plot H3K27me3 coverage fold change mean or median (+Dox/Ctrl) [y] over fold expression change mean or median (+Dox/Ctrl) log[x].
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    K27exprCorrDF=pd.DataFrame(columns=["chrom","start","WT_expr_averageFold", "WT_TSS_averageFold","WT_Gene_averageFold","KO_expr_averageFold", "KO_TSS_averageFold","KO_Gene_averageFold"])
    AverageFoldDF=getAverageFold(analysisDF, calculate_median=calculate_median)       # get average fold expression change on Xist induction from WT and Ubn1/2 KO
    K27exprCorrDF[["chrom", "start", "WT_expr_averageFold", "KO_expr_averageFold"]]=AverageFoldDF[["chrom","start","WT","Ubn1/2"]]

    cols=[]
    for gtp in ["WT", "KO"]:
        for repl in ["N8","N9","N10"]:
            cols.append(gtp+"_"+repl+TSS_suffix)
            cols.append(gtp+"_"+repl+Gene_suffix)
    ChIPfoldDF=pd.DataFrame(columns=cols)
    for gtp in ["WT", "KO"]:
        for repl in ["N8","N9","N10"]:
            ChIPfoldDF[gtp+"_"+repl+TSS_suffix]=ChIP_DF[gtp+"_DOX_"+repl+TSS_suffix]/ChIP_DF[gtp+"_"+repl+TSS_suffix]
            ChIPfoldDF[gtp+"_"+repl+Gene_suffix]=ChIP_DF[gtp+"_DOX_"+repl+Gene_suffix]/ChIP_DF[gtp+"_"+repl+Gene_suffix]        
    for gtp in ["WT", "KO"]:
        TSS_cols=[]
        Gene_cols=[]
        for repl in ["N8","N9","N10"]:
            TSS_cols.append(gtp+"_"+repl+TSS_suffix)
            Gene_cols.append(gtp+"_"+repl+Gene_suffix)
        if calculate_median:
            K27exprCorrDF[gtp+"_TSS_averageFold"]=ChIPfoldDF[TSS_cols].median(axis=1)
            K27exprCorrDF[gtp+"_Gene_averageFold"]=ChIPfoldDF[Gene_cols].median(axis=1)
        else:
            K27exprCorrDF[gtp+"_TSS_averageFold"]=ChIPfoldDF[TSS_cols].mean(axis=1)
            K27exprCorrDF[gtp+"_Gene_averageFold"]=ChIPfoldDF[Gene_cols].mean(axis=1)
    K27exprCorrDF=K27exprCorrDF.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    K27exprCorrDF=K27exprCorrDF.dropna(axis=0, how="any")
    X_K27expFoldDF=K27exprCorrDF[ K27exprCorrDF["chrom"]=="chrX" ]
    A_K27expFoldDF=K27exprCorrDF[ K27exprCorrDF["chrom"]!="chrX" ]
    fig=plt.figure(num="K27_expr_fold", figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])
    panel_row=0
    for gtp in ["WT","KO"]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 0])
            ax.set_title(gtp+" fold expression vs K27 TSS")
            ax.set_xscale("log")
            ax.plot(A_K27expFoldDF[gtp+"_expr_averageFold"].values, A_K27expFoldDF[gtp+"_TSS_averageFold"].values, ".", color="darkcyan")
            ax.plot(X_K27expFoldDF[gtp+"_expr_averageFold"].values, X_K27expFoldDF[gtp+"_TSS_averageFold"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
            ax=fig.add_subplot(gs[panel_row, 1])
            ax.set_title(gtp+" fold expression vs K27 gene body")
            ax.set_xscale("log")
            ax.plot(A_K27expFoldDF[gtp+"_expr_averageFold"].values, A_K27expFoldDF[gtp+"_Gene_averageFold"].values, ".", color="darkcyan")
            ax.plot(X_K27expFoldDF[gtp+"_expr_averageFold"].values, X_K27expFoldDF[gtp+"_Gene_averageFold"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
        panel_row+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("K27_expr_fold")

def K27foldExpressionLevelPlot(analysis, ChIP_DF, calculate_median=False, filename="K27_fold_exprLevel.png"):
    """
    Plot H3K27me3 coverage fold change mean or median (+Dox/Ctrl) [y] over mean or median expression level in Ctrl samples log[x].
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    K27exprCorrDF=pd.DataFrame(columns=["chrom","start","expression_level", "WT_TSS_averageFold","WT_Gene_averageFold", "KO_TSS_averageFold","KO_Gene_averageFold"])
    K27exprCorrDF["expression_level"]=getAverageExpressionLevel(analysisDF, calculate_median=calculate_median)       # get average fold expression change on Xist induction from WT and Ubn1/2 KO
    K27exprCorrDF[["chrom", "start"]]=analysisDF[["chrom","start"]]

    cols=[]
    for gtp in ["WT", "KO"]:
        for repl in ["N8","N9","N10"]:
            cols.append(gtp+"_"+repl+TSS_suffix)
            cols.append(gtp+"_"+repl+Gene_suffix)
    ChIPfoldDF=pd.DataFrame(columns=cols)
    for gtp in ["WT", "KO"]:
        for repl in ["N8","N9","N10"]:
            ChIPfoldDF[gtp+"_"+repl+TSS_suffix]=ChIP_DF[gtp+"_DOX_"+repl+TSS_suffix]/ChIP_DF[gtp+"_"+repl+TSS_suffix]
            ChIPfoldDF[gtp+"_"+repl+Gene_suffix]=ChIP_DF[gtp+"_DOX_"+repl+Gene_suffix]/ChIP_DF[gtp+"_"+repl+Gene_suffix]        
    for gtp in ["WT", "KO"]:
        TSS_cols=[]
        Gene_cols=[]
        for repl in ["N8","N9","N10"]:
            TSS_cols.append(gtp+"_"+repl+TSS_suffix)
            Gene_cols.append(gtp+"_"+repl+Gene_suffix)
        if calculate_median:
            K27exprCorrDF[gtp+"_TSS_averageFold"]=ChIPfoldDF[TSS_cols].median(axis=1)
            K27exprCorrDF[gtp+"_Gene_averageFold"]=ChIPfoldDF[Gene_cols].median(axis=1)
        else:
            K27exprCorrDF[gtp+"_TSS_averageFold"]=ChIPfoldDF[TSS_cols].mean(axis=1)
            K27exprCorrDF[gtp+"_Gene_averageFold"]=ChIPfoldDF[Gene_cols].mean(axis=1)
    K27exprCorrDF=K27exprCorrDF.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    K27exprCorrDF=K27exprCorrDF.dropna(axis=0, how="any")
    X_K27expFoldDF=K27exprCorrDF[ K27exprCorrDF["chrom"]=="chrX" ]
    A_K27expFoldDF=K27exprCorrDF[ K27exprCorrDF["chrom"]!="chrX" ]
    fig=plt.figure(num="K27_fold_exprLevel", figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])
    panel_row=0
    for gtp in ["WT","KO"]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 0])
            ax.set_title(gtp+" TSS K27 fold vs expression level")
            ax.set_xscale("log")
            ax.plot(A_K27expFoldDF["expression_level"].values, A_K27expFoldDF[gtp+"_TSS_averageFold"].values, ".", color="darkcyan")
            ax.plot(X_K27expFoldDF["expression_level"].values, X_K27expFoldDF[gtp+"_TSS_averageFold"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
            ax=fig.add_subplot(gs[panel_row, 1])
            ax.set_title(gtp+" gene body K27 fold vs expression level")
            ax.set_xscale("log")
            ax.plot(A_K27expFoldDF["expression_level"].values, A_K27expFoldDF[gtp+"_Gene_averageFold"].values, ".", color="darkcyan")
            ax.plot(X_K27expFoldDF["expression_level"].values, X_K27expFoldDF[gtp+"_Gene_averageFold"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
        panel_row+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("K27_fold_exprLevel")

def XposK27FoldPlot(analysis, ChIP_DF, calculate_median=False, filename="Xpos_K27_fold.png"):
    """
    Plot H3K27me3 coverage fold change mean or median (+Dox/Ctrl) [y] over genomic position on the X chromosome.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    K27exprCorrDF=pd.DataFrame(columns=["chrom","start","WT_TSS_averageFold","WT_Gene_averageFold","KO_TSS_averageFold","KO_Gene_averageFold"])
    K27exprCorrDF[["chrom", "start"]]=analysisDF[["chrom","start"]]
    cols=[]
    for gtp in ["WT", "KO"]:
        for repl in ["N8","N9","N10"]:
            cols.append(gtp+"_"+repl+TSS_suffix)
            cols.append(gtp+"_"+repl+Gene_suffix)
    ChIPfoldDF=pd.DataFrame(columns=cols)
    for gtp in ["WT", "KO"]:
        for repl in ["N8","N9","N10"]:
            ChIPfoldDF[gtp+"_"+repl+TSS_suffix]=ChIP_DF[gtp+"_DOX_"+repl+TSS_suffix]/ChIP_DF[gtp+"_"+repl+TSS_suffix]
            ChIPfoldDF[gtp+"_"+repl+Gene_suffix]=ChIP_DF[gtp+"_DOX_"+repl+Gene_suffix]/ChIP_DF[gtp+"_"+repl+Gene_suffix]        
    for gtp in ["WT", "KO"]:
        TSS_cols=[]
        Gene_cols=[]
        for repl in ["N8","N9","N10"]:
            TSS_cols.append(gtp+"_"+repl+TSS_suffix)
            Gene_cols.append(gtp+"_"+repl+Gene_suffix)
        if calculate_median:
            K27exprCorrDF[gtp+"_TSS_averageFold"]=ChIPfoldDF[TSS_cols].median(axis=1)
            K27exprCorrDF[gtp+"_Gene_averageFold"]=ChIPfoldDF[Gene_cols].median(axis=1)
        else:
            K27exprCorrDF[gtp+"_TSS_averageFold"]=ChIPfoldDF[TSS_cols].mean(axis=1)
            K27exprCorrDF[gtp+"_Gene_averageFold"]=ChIPfoldDF[Gene_cols].mean(axis=1)
    K27exprCorrDF=K27exprCorrDF.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    K27exprCorrDF=K27exprCorrDF.dropna(axis=0, how="any")
    X_K27expFoldDF=K27exprCorrDF[ K27exprCorrDF["chrom"]=="chrX" ]
    fig=plt.figure(num="Xpos_K27_fold", figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])
    xticks=list(range(0,max(X_K27expFoldDF["start"].values),20000000))
    xtick_labels=list(range(0,max(xticks)/1000000+1,20))
    panel_row=0
    for gtp in ["WT","KO"]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 0])         # plot TSS K27 coverage
            ax.set_title(gtp+" TSS fold K27 over chrX pos")
            ax.plot(X_K27expFoldDF["start"].values, X_K27expFoldDF[gtp+"_TSS_averageFold"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels)
            ax.tick_params(which='major', length=4, color='black')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 1])         # plot gene body K27 coverage
            ax.set_title(gtp+" gene body fold K27 over chrX pos")
            ax.plot(X_K27expFoldDF["start"].values, X_K27expFoldDF[gtp+"_Gene_averageFold"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels)
            ax.tick_params(which='major', length=4, color='black')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        panel_row+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("Xpos_K27_fold")

def K27exprCorrelationPlot(analysis, ChIP_DF, calculate_median=False, filename="K27_expr_correlation.png"):
    """
    Plot H3K27me3 coverage [y] over expression level of the genes in the Ctrl samples.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    K27exprCorrDF=pd.DataFrame(columns=["chrom","start","expression_level", "WT_TSS_K27avg", "WT_Gene_K27avg", "KO_TSS_K27avg","KO_Gene_avg"])
    K27exprCorrDF["expression_level"]=getAverageExpressionLevel(analysisDF, calculate_median=calculate_median)       # get average fold expression change on Xist induction from WT and Ubn1/2 KO
    K27exprCorrDF[["chrom", "start"]]=analysisDF[["chrom","start"]]
    for gtp in ["WT", "KO"]:
        TSS_cols=[]
        Gene_cols=[]
        for repl in ["N8","N9","N10"]:
            TSS_cols.append(gtp+"_"+repl+TSS_suffix)
            Gene_cols.append(gtp+"_"+repl+Gene_suffix)
        if calculate_median:
            K27exprCorrDF[gtp+"_TSS_K27avg"]=ChIP_DF[TSS_cols].median(axis=1)
            K27exprCorrDF[gtp+"_Gene_K27avg"]=ChIP_DF[Gene_cols].median(axis=1)
        else:
            K27exprCorrDF[gtp+"_TSS_K27avg"]=ChIP_DF[TSS_cols].mean(axis=1)
            K27exprCorrDF[gtp+"_Gene_K27avg"]=ChIP_DF[Gene_cols].mean(axis=1)
    X_K27expDF=K27exprCorrDF[ K27exprCorrDF["chrom"]=="chrX" ]
    A_K27expDF=K27exprCorrDF[ K27exprCorrDF["chrom"]!="chrX" ]
    fig=plt.figure(num="K27exprCorr", figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,1], width_ratios=[1,1])
    panel_row=0
    for gtp in ["WT","KO"]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, 0])
            ax.set_title(gtp+" TSS K27 expression correlation")
            ax.set_xscale("log")
            ax.plot(A_K27expDF["expression_level"].values, A_K27expDF[gtp+"_TSS_K27avg"].values, ".", color="darkcyan")
            ax.plot(X_K27expDF["expression_level"].values, X_K27expDF[gtp+"_TSS_K27avg"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
            ax=fig.add_subplot(gs[panel_row, 1])
            ax.set_title(gtp+" gene body K27 expression correlation")
            ax.set_xscale("log")
            ax.plot(A_K27expDF["expression_level"].values, A_K27expDF[gtp+"_Gene_K27avg"].values, ".", color="darkcyan")
            ax.plot(X_K27expDF["expression_level"].values, X_K27expDF[gtp+"_Gene_K27avg"].values, ".", color="firebrick")
            ax.set_ylim(bottom=0, top=50)
        panel_row+=1
    fig.tight_layout()
    fig.savefig(DATAPATH+filename)
    plt.close("K27exprCorr")

def HeizMappe(analysis, column_order=None, LinReg=True, chrom="chr1", colorMap=None, genotype_sepLine_color=None, filename=None):
    """
    Generated a pyplot figure using seaborn for X-chromosomal gene expression and a reference autosome. If a filename is given the figure is saved else displayed
    using matplotlib. It accepts an analysis DataFrame or a XpressionAnalysis object. Normalization for the X chromosome can be done by calculating
    a parameter for each cell line in Ctrl conditions and applying it to the cell lines in Ctrl and Dox conditions. This will rrequire that the global variables
    genotypes, conditions, and samples are set.

    Parameters:

        LinReg (boolean)      : Specified to use linear regression to adjust for X-linked gene expression between the samples using the no dox samples for estimating
                               the number of X-chromosomes for both dox and no dox datasets of one sample pair
        chrom (str)           : The reference chromosome to display as a control (default = chr1)
        filename (str, None)  : A valid file path to save the figure to. If None the figure will be displayed
        column_order          : A list of samples in the order of the columns of the resulting heatmap. If None a default sample order will be generated.
        colorMap              : Specifies the color map for the heatmap
        genotype_sepLine_color: If not None separating lines will be drawn with this color in between genotypes assuming 2 conditions and 3 replicates per genotype

    The function presumes that samples, genotypes, and conditions ("Ctrl", "Dox") are initialized as global variables.        
    """
    if colorMap is None:
        colorMap=XpressionAnalysis.RBGcmap
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    XgeneDF=analysisDF[ analysisDF["chrom"]=="chrX" ]
    AgeneDF=analysisDF[ analysisDF["chrom"]==chrom ]
    XgeneDF=XgeneDF.sort_values( by=["start"], ascending=True)
    AgeneDF=AgeneDF.sort_values( by=["start"], ascending=True)
    XgeneDFnorm=pd.DataFrame(columns=XgeneDF.columns) 
    XgeneDFnorm[XpressionAnalysis.gene_annotation]=XgeneDF[XpressionAnalysis.gene_annotation]    # copy the gene annotation
    slopes={}
    ref_replicate=None
    for genotype in genotypes:
        for replicate in samples[genotype]["Ctrl"]:
            if ref_replicate is None:
                slopes[replicate]=1                                  # set the reference sample slope to 1
                ref_replicate=replicate
            else:
                slopes[replicate]=XpressionAnalysis.lin_reg_slope_no_intercept(XgeneDF[ref_replicate].values, XgeneDF[replicate].values) # get slope of linear regression line between replicate and reference
    for genotype in genotypes:
        for replIdx in range(len(samples[genotype]["Ctrl"])):
            replCtrl=samples[genotype]["Ctrl"][replIdx]
            replDox=samples[genotype]["Dox"][replIdx]
            XgeneDFnorm[replCtrl]=XgeneDF[replCtrl]/slopes[replCtrl]
            XgeneDFnorm[replDox]=XgeneDF[replDox]/slopes[replCtrl]    # normalize both the Ctrl and Dox data with the slope of the Ctrl sample of the replicate
    if column_order is None:
        column_order=[]
        for genotype in ["WT","Ubn1","Ubn2","Hira","Ubn1/2"]:         # genotypes:
            for condition in conditions:
                for replicate in samples[genotype][condition]:
                    column_order.append(replicate)
    XgeneDFnorm.to_csv(DATAPATH+"chrX_expressionDF_normalizedToCtrl.csv")
    XgeneDF.to_csv(DATAPATH+"chrX_expressionDF.csv")
    AgeneDF.to_csv(DATAPATH+chrom+"_expressionDF.csv")
    exprX=XgeneDF[column_order].values
    exprA=AgeneDF[column_order].values
    exprXnorm=XgeneDFnorm[column_order].values
    zmapX=scipy.stats.zscore(exprX, axis=1)
    zmapX=zmapX[~np.isnan(zmapX).any(axis=1)]
    zmapA=scipy.stats.zscore(exprA, axis=1)
    zmapA=zmapA[~np.isnan(zmapA).any(axis=1)]
    zmapXnorm=scipy.stats.zscore(exprXnorm, axis=1)
    zmapXnorm=zmapXnorm[~np.isnan(zmapXnorm).any(axis=1)]
    fig=plt.figure(num="RNAseq_Heatmap", facecolor='w', edgecolor='k')   # , figsize=(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1, ncols=4, width_ratios=[18,18,18,1])
    ax1=fig.add_subplot(gs[0, 0])
    ax2=fig.add_subplot(gs[0, 1])
    ax3=fig.add_subplot(gs[0, 2])
    ax4=fig.add_subplot(gs[0, 3])
    ax1.set_title("chrX", loc='left', y=1.08)
    ax1.title.set_fontsize(20)
    seaborn.heatmap(zmapX, cmap=colorMap, xticklabels=False, yticklabels=False, ax=ax1, cbar=False)
    ax1.xaxis.set_tick_params(labeltop='on', labelbottom='off')
    ax1.set_xticklabels(column_order)
    for label in ax1.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(12)
    if not (genotype_sepLine_color is None):
        for x in range(len(column_order)/3-2):
            ax1.axvline(x=(x+1)*6, linewidth=1, color=genotype_sepLine_color) 
    ax2.set_title("chr1", loc='left', y=1.08)
    ax2.title.set_fontsize(20)
    seaborn.heatmap(zmapA, cmap=colorMap, xticklabels=False, yticklabels=False, ax=ax2, cbar=False)
    ax2.xaxis.set_tick_params(labeltop='on', labelbottom='off')
    ax2.set_xticklabels(column_order)
    for label in ax2.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(12)
    if not (genotype_sepLine_color is None):
        for x in range(len(column_order)/3-2):
            ax2.axvline(x=(x+1)*6, linewidth=1, color=genotype_sepLine_color) 
    ax3.set_title("chrX normalized Xs", loc='left', y=1.08)
    ax3.title.set_fontsize(20)
    seaborn.heatmap(zmapXnorm, cmap=colorMap, xticklabels=False, yticklabels=False, ax=ax3, cbar=True, cbar_ax=ax4)
    ax3.xaxis.set_tick_params(labeltop='on', labelbottom='off')
    ax3.set_xticklabels(column_order)
    for label in ax3.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(12)
    if not (genotype_sepLine_color is None):
        for x in range(len(column_order)/3-2):
            ax3.axvline(x=(x+1)*6, linewidth=1, color=genotype_sepLine_color) 
    if (filename is None):
        plt.show()
    else:
        fig.savefig(DATAPATH+filename)
    plt.close("RNAseq_Heatmap")

def plotGenes(analysis, genes, genotypes, filename=None):
    """
    Plots the read or TPM counts for a list of genes for the genotypes indicated.

    analysis          : XpresssionAnalysis object or DataFrame with the RNAseq analysis dataset
    genes (list)      : Gene names for which expression is plotted
    filename (string) : Name of file to which the resulting figure is saved. If None the figure will be shown on screen

    The function presumes that samples, and conditions are initialized as global variables.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    color={"Ctrl": "black", "Dox": "firebrick"}
    fig=plt.figure(num="plotGene", facecolor='w', edgecolor='k')   # , figsize=(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1+int((len(genes)-0.5)/3), ncols=3, width_ratios=[1]*3)
    panel_row=0
    panel_col=0
    for gene in genes:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(gene)
            xlabs=[]
            x=0.5
            for genotype in genotypes:
                xlabs.append(genotype)
                for condition in conditions:
                    v=[]
                    for replicate in samples[genotype][condition]:
                        v.append(analysisDF.loc[gene,replicate])
                        ax.plot(x, analysisDF.loc[gene,replicate], "o", color=color[condition])
                    ax.hlines(np.mean(v),x-0.3,x+0.3, color=color[condition])                        # mean line
                    ax.hlines(np.mean(v)-np.std(v),x-0.15,x+0.15, color=color[condition])            # -std line
                    ax.hlines(np.mean(v)+np.std(v),x-0.15,x+0.15, color=color[condition])            # -std line
                    ax.vlines(x,np.mean(v)-np.std(v),np.mean(v)+np.std(v), color=color[condition])   # vertical line -std to +std
                    x+=1
            ax.set_ylim(bottom=0)
            xticks=[2*i+1 for i in range(len(genotypes))]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabs)
        panel_col+=1
        if panel_col>=3:
            panel_col=0
            panel_row+=1
    fig.tight_layout()
    if (filename is None):
        plt.show()
    else:
        fig.savefig(DATAPATH+filename)    
    plt.close("plotGene")

def plotGenesFoldChange(analysis, genes, genotypes, filename=None, colorWT="green"):
    """
    Plots the fold change between "Dox" and "Ctrl" from read or TPM counts for a list of genes for the genotypes indicated.

    analysis          : XpresssionAnalysis object or DataFrame with the RNAseq analysis dataset
    genes (list)      : Gene names for which expression is plotted
    filename (string) : Name of file to which the resulting figure is saved. If None the figure will be shown on screen

    The function presumes that samples, and conditions are initialized as global variables.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    color={"Ctrl": "black", "Dox": "firebrick"}
    fig=plt.figure(num="plotGenesFoldChange", facecolor='w', edgecolor='k')   # , figsize=(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1+int((len(genes)-0.5)/3), ncols=3, width_ratios=[1]*3)
    panel_row=0
    panel_col=0
    for gene in genes:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(gene)
            xlabs=[]
            x=0.5
            for genotype in genotypes:
                v={ "Ctrl": [], "Dox": [] }                       # initialize a dictionary for collecting read or TPM counts for each sample and condition
                for condition in conditions:
                    for replicate in samples[genotype][condition]:
                        v[condition].append(analysisDF.loc[gene,replicate])
                fold_change = []
                for d, c in zip(v["Dox"], v["Ctrl"]):
                    if c!=0:                                     # eliminate values when control has 0 counts avoiding div/0
                        fold_change.append(d/c)
                if len(fold_change)>0:
                    if genotype=="WT":
                        color=colorWT
                        ax.axhline(y=np.mean(fold_change), ls="--", color=colorWT, alpha=0.5)
                    else:
                        color="black"
                    ax.plot([x]*len(fold_change), fold_change, "o", color=color)
                    ax.hlines(np.mean(fold_change),x-0.3,x+0.3, color=color)                                                      # mean line
                    ax.hlines(np.mean(fold_change)-np.std(fold_change),x-0.15,x+0.15, color=color)                                # -std line
                    ax.hlines(np.mean(fold_change)+np.std(fold_change),x-0.15,x+0.15, color=color)                                # +std line
                    ax.vlines(x,np.mean(fold_change)-np.std(fold_change),np.mean(fold_change)+np.std(fold_change), color=color)   # vertical line -std to +std
                    xlabs.append(genotype)
                    x+=1
            y0, y1 = ax.get_ylim()
            y1 = max(1.2,y1)
            ax.set_ylim(bottom=0, top=y1)
            xticks=[i+0.5 for i in range(len(genotypes))]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabs)
            ax.axhline(y=1.0, linestyle="--", color="black", alpha=0.5)
        panel_col+=1
        if panel_col>=3:
            panel_col=0
            panel_row+=1
    fig.tight_layout()
    if (filename is None):
        plt.show()
    else:
        fig.savefig(DATAPATH+filename)    
    plt.close("plotGenesFoldChange")

def differentialGeneExpressionAnalysis(analysis, refGrp, sbjGrp, fold_change_threshold=2.0, p_value_cutoff=0.01, chrom=None, normalize=None, filename=None, ax=None, show_figure=True, show_Xist=False):
    """
    Plots the fold change between "Dox" and "Ctrl" from read or TPM counts for a list of genes for the genotypes indicated.

    analysis          : XpresssionAnalysis object or DataFrame with the RNAseq analysis dataset
    genes (list)      : Gene names for which expression is plotted
    filename (string) : Name of file to which the resulting figure is saved. If None the figure will be shown on screen

    The function presumes that samples, and conditions are initialized as global variables.
    """
    if type(analysis) == pd.DataFrame:
        analysisDF=analysis
    else:
        analysisDF=analysis.getDataFrame()
    if not (chrom is None):
        analysisDF=analysisDF[ analysisDF["chrom"]==chrom ]
    analysisDF["ref_mean"]=analysisDF[refGrp].mean(axis=1)
    analysisDF["sbj_mean"]=analysisDF[sbjGrp].mean(axis=1)
    analysisDF["fold"]=analysisDF["sbj_mean"]/analysisDF["ref_mean"]
    analysisDF["p-value"]=scipy.stats.ttest_ind(analysisDF[refGrp], analysisDF[sbjGrp], axis=1).pvalue # [1]
    # analysisDF.to_csv(DATAPATH+"DifExpr_"+str(chrom)+".csv")
    fig=None
    if ax is None:
        fig=plt.figure(num="differentialGeneExpression", facecolor='w', edgecolor='k')
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(111)
    if not(filename is None):
        ax.set_title(filename.split(".")[0])
    ax.set_xscale("log")
    ax.set_yscale("log")
    filteredDF=analysisDF.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    filteredDF=filteredDF.dropna(axis=0, how="any")
    ax.plot(filteredDF["fold"], 1/filteredDF["p-value"], ".", color="gray", alpha=0.5)
    analysisDF=analysisDF.sort_values(by="fold", ascending=False)
    upDF=analysisDF[ analysisDF["fold"]>=fold_change_threshold ]    
    upDF=upDF[ upDF["p-value"]<=p_value_cutoff ]
    downDF=analysisDF[ analysisDF["fold"]<=1.0/float(fold_change_threshold) ]    
    downDF=downDF[ downDF["p-value"]<=p_value_cutoff ]
    if not (filename is None):
        analysisDF.to_csv(DATAPATH+filename+".csv")
    upGenesList=list(upDF.index)
    downGenesList=list(downDF.index)
    filteredDF=upDF.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    filteredDF=filteredDF.dropna(axis=0, how="any")
    ax.plot(filteredDF["fold"], 1/filteredDF["p-value"], ".", color="darkcyan")
    if show_Xist and ("Xist" in filteredDF.index):
        ax.plot(filteredDF.loc["Xist"]["fold"], 1/filteredDF.loc["Xist"]["p-value"], "o", color="darkcyan")
        ax.text(filteredDF.loc["Xist"]["fold"]*0.7, 2/filteredDF.loc["Xist"]["p-value"], "Xist", color="black")
    filteredDF=downDF.replace([np.inf, -np.inf], np.nan)   # filter infs and NaNs
    filteredDF=filteredDF.dropna(axis=0, how="any")
    ax.plot(filteredDF["fold"], 1/filteredDF["p-value"], ".", color="firebrick")
    if not (fig is None):
        fig.tight_layout()
    if show_figure:
        plt.show()
    if not (filename is None):
        upDF.to_csv(DATAPATH+filename.split(".")[0]+"UP.csv")
        upDF.to_csv(DATAPATH+filename.split(".")[0]+"DOWN.csv")
        fig.savefig(DATAPATH+filename)    
    if not (fig is None):
        plt.close("differentialGeneExpression")
    print len(upGenesList),len(upDF),"up and", len(downGenesList), len(downDF), "down regulated"
    return upGenesList, downGenesList

def difGeneExprAllGenotypes(analysis, genotypes, filename=None, fold_change_threshold=2.0, p_value_cutoff=0.01):
    difRegGenes={}
    if len(genotypes)<2:
        if printVerboseLog:
            print "ERROR difGeneExprAllGenotypes: requires at least 2 genotypes to compare.", len(genotypes), "provided."
        return None
    fig=plt.figure(num="difGenesExpr", facecolor='w', edgecolor='k')   # , figsize=(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1+int((len(genotypes)-1.5)/2), ncols=min(len(genotypes)-1,2), width_ratios=[1]*2, height_ratios=[1]*(1+int((len(genotypes)-1.5)/2)))
    panel_row=0
    panel_col=0
    refGt=genotypes[0]
    for gt in genotypes[1:]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(gt+" vs. "+refGt)
            difRegGenes[gt]=dict()
            difRegGenes[gt]["up"], difRegGenes[gt]["down"]=differentialGeneExpressionAnalysis(analysis, samples[refGt]["Ctrl"], samples[gt]["Ctrl"], filename=None, fold_change_threshold=fold_change_threshold, p_value_cutoff=p_value_cutoff, ax=ax, show_figure=False)
        panel_col+=1
        if panel_col>=2:
            panel_col=0
            panel_row+=1
    fig.tight_layout()
    if (filename is None):
        plt.show()
    else:
        fig.savefig(DATAPATH+filename)    
    plt.close("difGenesExpr")
    return difRegGenes

def Dox_vs_Ctrl_AllGenotypes(analysis, genotypes, filename=None, show_Xist=False, fold_change_threshold=2.0, p_value_cutoff=0.01):
    fig=plt.figure(num="Dox_vs_Ctrl", facecolor='w', edgecolor='k')   # , figsize=(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1+int((len(genotypes)-0.5)/2), ncols=min(len(genotypes)-1,2), width_ratios=[1]*2, height_ratios=[1]*(1+int((len(genotypes)-0.5)/2)))
    panel_row=0
    panel_col=0
    for gt in genotypes:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_title(gt)
            differentialGeneExpressionAnalysis(analysis, samples[gt]["Ctrl"], samples[gt]["Dox"], filename=None, ax=ax, fold_change_threshold=fold_change_threshold, p_value_cutoff=p_value_cutoff, show_figure=False, show_Xist=show_Xist)
        panel_col+=1
        if panel_col>=2:
            panel_col=0
            panel_row+=1
    fig.tight_layout()
    if (filename is None):
        plt.show()
    else:
        fig.savefig(DATAPATH+filename)    
    plt.close("Dox_vs_Ctrl")

def plotVennDiagrams(difRegGenes, genotypes, filename=None):
    fig=plt.figure(num="VennDiagrams", figsize=(10, 6), dpi=300, facecolor='w', edgecolor='k')   # , figsize=(20, 12), dpi=300
    gs=gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1]*2)
    panel_col=0
    for reg in ["up", "down"]:
        with plt.style.context('seaborn-white'):
            ax=fig.add_subplot(gs[0, panel_col])
            ax.set_title(reg+"regulated genes")
            sets=[]
            names=[]
            for gt in genotypes:
                sets.append(set(difRegGenes[gt][reg]))
                names.append(gt)                            
            venn3(sets, names)
        panel_col+=1
    fig.tight_layout()
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)
    plt.close("VennDiagrams")

# K27ExpFoldPlot(tpmDF, ChIP_DF)
# XposK27FoldPlot(tpmDF, ChIP_DF)
# K27foldExpressionLevelPlot(tpmDF, ChIP_DF)
# fexGroupBarPlot(tpmDF, calculate_median=False, filename="fold_expression_correlation.pdf")
# K27exprCorrelationPlot(tpmDF, ChIP_DF, calculate_median=False, filename="K27_expr_correlation.png")
# averageFoldHistogram(foldFilteredDF, bins=20, smooth=False, density=True, show_median_lines=True, filename="XA_averageFold_histogram.pdf")
# averageFoldHistogram(foldFilteredDF, bins=20, smooth=True, density=True, show_median_lines=True, filename="smooth_XA_averageFold_histogram.pdf")
# tpmAnalysis=XpressionAnalysis.from_df(tpmDF)
# tpmAnalysis.HeizMappe(getPairedSampleList())
# HeizMappe(tpmDF, colorMap="seismic", genotype_sepLine_color="#000000")
# plotGenes(tpmDF, ["Xist","Pgk1","Hprt","Mecp2","Rlim","Atrx","Mid1","Kdm6a","Kdm5c","Pou5f1","Sox2","Gapdh"],["WT","Ubn1","Ubn2","Hira","Ubn1/2"])
# plotGenes(tpmDF, ["Ubn1","Ubn2","Hira","Pou5f1","Sox2","Gapdh"],["WT","Ubn1","Ubn2","Hira","Ubn1/2"])
# plotGenesFoldChange(tpmDF, ["Xist","Pgk1","Hprt","Mecp2","Rlim","Atrx","Mid1","Kdm6a","Kdm5c","Pou5f1","Sox2","Gapdh"],["WT","Ubn1/2"])
# plotGenesFoldChange(tpmDF, ["Xist","Pgk1","Hprt","Mecp2","Rlim","Atrx","Mid1","Kdm6a","Kdm5c","Pou5f1","Sox2","Gapdh"],["WT","Ubn1","Ubn2","Hira","Ubn1/2"])
# differentialGeneExpressionAnalysis(tpmDF, samples["WT"]["Ctrl"], samples["Ubn1/2"]["Ctrl"])
# differentialGeneExpressionAnalysis(tpmDF, samples["WT"]["Ctrl"], samples["WT"]["Dox"])
# differentialGeneExpressionAnalysis(tpmDF, samples["Ubn1/2"]["Ctrl"], samples["Ubn1/2"]["Dox"])
# difGeneExprAllGenotypes(tpmDF, ["WT","Ubn1","Ubn2","Hira","Ubn1/2"])
# Dox_vs_Ctrl_AllGenotypes(tpmDF, ["WT","Ubn1","Ubn2","Hira","Ubn1/2"], show_Xist=True)

# difRegGenes=difGeneExprAllGenotypes(tpmDF, ["WT","Ubn1","Ubn2","Hira","Ubn1/2"], fold_change_threshold=2, p_value_cutoff=0.01)
# plotVennDiagrams(difRegGenes,["Ubn1","Ubn2","Ubn1/2"], filename="VennDiagram_Ubn1_Ubn2_Ubn12.pdf")
# plotVennDiagrams(difRegGenes,["Ubn2","Hira","Ubn1/2"], filename="VennDiagram_Ubn2_Hira_Ubn12.pdf")
"""
for gt in difRegGenes.keys():
    with open("UpregulatedGenes_"+gt.replace("/","_")+".txt","w") as outfile:
        for gene in difRegGenes[gt]["up"]:
            outfile.writelines(gene+"\n")
    with open("DownregulatedGenes_"+gt.replace("/","_")+".txt","w") as outfile:
        for gene in difRegGenes[gt]["down"]:
            outfile.writelines(gene+"\n")
"""
















