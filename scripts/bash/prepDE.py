#!/usr/bin/env python3
import re, csv, sys, os, glob, warnings, itertools
from math import ceil
from optparse import OptionParser
from operator import itemgetter
from collections import defaultdict

parser=OptionParser(description='Generates two CSV files containing the count matrices for genes and transcripts, using the coverage values found in the output of `stringtie -e` or TACO')
parser.add_option('-i', '--input', '--in', default='.', help="a folder containing all sample sub-directories, or a text file with sample ID and path to its GTF file on each line [default: %default/]")
parser.add_option('-g', default='gene_count_matrix.csv', help="where to output the gene count matrix [default: %default")
parser.add_option('-t', default='transcript_count_matrix.csv', help="where to output the transcript count matrix [default: %default]")
parser.add_option('-l', '--length', default=75, type='int', help="the average read length [default: %default]")
parser.add_option('-p', '--pattern', default=".", help="a regular expression that selects the sample subdirectories")
parser.add_option('-c', '--cluster', action="store_true", help="whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)")
parser.add_option('-s', '--string', default="MSTRG", help="if a different prefix is used for geneIDs assigned by StringTie [default: %default]")
parser.add_option('-k', '--key', default="prepG", help="if clustering, what prefix to use for geneIDs assigned by this script [default: %default]")
parser.add_option('-v', action="store_true", help="enable verbose processing")

parser.add_option('--legend', default="legend.csv", help="if clustering, where to output the legend file mapping transcripts to assigned geneIDs [default: %default]")
(opts, args)=parser.parse_args()

samples = [] # List of tuples. If sample list, (first column, path). Else, (subdirectory name, path to gtf file in subdirectory)
if (os.path.isfile(opts.input)):
    # gtfList = True
    try:
        fin = open(opts.input, 'r')
        for line in fin:
            if line[0] != '#':
                lineLst = tuple(line.strip().split(None,2))
                if (len(lineLst) != 2):
                    print("Error: line should have a sample ID and a file path:\n%s" % (line.strip()))
                    exit(1)
                if lineLst[0] in samples:
                    print("Error: non-unique sample ID (%s)" % (lineLst[0]))
                    exit(1)
                if not os.path.isfile(lineLst[1]):
                    print("Error: GTF file not found (%s)" % (lineLst[1]))
                    exit(1)
                samples.append(lineLst)
    except IOError:
        print("Error: List of .gtf files, %s, doesn't exist" % (opts.input))
        exit(1)
else:
    # gtfList = False
    ## Check that opts.input directory exists
    if not os.path.isdir(opts.input):
      parser.print_help()
      print(" ")
      print("Error: sub-directory '%s' not found!" % (opts.input))
      sys.exit(1)
    #####
    ## Collect all samples file paths and if empty print help message and quit
    #####
    samples = []
    for i in next(os.walk(opts.input))[1]:
        if re.search(opts.pattern,i):
         for f in glob.iglob(os.path.join(opts.input,i,"*.gtf")):
            samples.append((i,f)) 

if len(samples) == 0:
  parser.print_help()
  print(" ")
  print("Error: no GTF files found under base directory %s !" % (opts.input))
  sys.exit(1)

# Regular expressions for parsing GTF attributes
# Modified to use raw strings (r'pattern') to avoid escape sequence warnings
RE_GENE_ID=re.compile(r'gene_id "([^"]+)"')
RE_GENE_NAME=re.compile(r'gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile(r'transcript_id "([^"]+)"')
RE_COVERAGE=re.compile(r'cov "([^"]+)"')
RE_STRING=re.compile(re.escape(opts.string))

# New regex patterns to handle TACO format
RE_TACO_TRANSCRIPT_ID=re.compile(r'ID=([^;]+)')
RE_TACO_GENE_ID=re.compile(r'geneID=([^;]+)')
RE_TACO_COVERAGE=re.compile(r'coverage=([^;\n]+)')

RE_GFILE=re.compile(r'\-G\s*(\S+)') # Use raw string to avoid escape warning


#####
## Sort the sample names by the sample ID
#####

samples.sort()

if opts.v:
  print("Sample GTFs found:")
  for s in samples:
     print("  %s: %s" % s)

#####
## Checks whether a given row is a transcript 
## other options: ex. exon, transcript, mRNA, 5'UTR
#####
def is_transcript(x):
  return len(x)>2 and x[2]=="transcript"

# Modified function to handle both formats of transcript IDs
def getTranscriptID(s):
    # Try StringTie format
    r=RE_TRANSCRIPT_ID.search(s)
    if r: 
        return r.group(1)
    # Try TACO format
    r=RE_TACO_TRANSCRIPT_ID.search(s)
    if r: 
        return r.group(1)
    return None

# Modified function to handle both formats of gene IDs
def getGeneID(s, ctg, tid):
    # Try StringTie format
    r=RE_GENE_ID.search(s)
    if r:
        rn=RE_GENE_NAME.search(s)
        if rn: 
            return r.group(1)+'|'+rn.group(1)
        else:
            return r.group(1)
    # Try TACO format
    r=RE_TACO_GENE_ID.search(s)
    if r: 
        return r.group(1)
    return tid

# Modified function to handle both formats of coverage values
def getCov(s):
    # Try StringTie format
    r=RE_COVERAGE.search(s)
    if r:
        try:
            v=float(r.group(1))
            if v<0.0: v=0.0
            return v
        except ValueError:
            pass
    # Try TACO format
    r=RE_TACO_COVERAGE.search(s)
    if r:
        try:
            v=float(r.group(1))
            if v<0.0: v=0.0
            return v
        except ValueError:
            pass
    return 0.0

def is_overlap(x,y): #NEEDS TO BE INTS!
  return x[0]<=y[1] and y[0]<=x[1]


def t_overlap(t1, t2): #from badGenes: chromosome, strand, cluster, start, end, (e1start, e1end)...
    if t1[0] != t2[0] or t1[1] != t2[1] or t1[5]<t2[4]: return False
    for i in range(6, len(t1)):
        for j in range(6, len(t2)):
            if is_overlap(t1[i], t2[j]): return True
    return False

## Average Readlength
read_len=opts.length

## Variables/Matrices to store t/g_counts
t_count_matrix, g_count_matrix=[],[]

##Get ready for clustering, stuff is once for all samples##
geneIDs=defaultdict(lambda: str) #key=transcript, value=cluster/gene_id


## For each of the sorted sample paths
for s in samples:
    badGenes=[] #list of bad genes (just ones that aren't MSTRG)
    try:
        ## opts.input = parent directory of sample subdirectories
        ## s = sample currently iterating through
        ## os.path.join(opts.input,s,"*.gtf") path to current sample's GTF
        ## split = list of lists: [[chromosome, ...],...]

        if opts.v:
            print("Processing file: %s" % s[1])
            
        # Read the GTF file line by line and filter out malformed lines
        split = []
        with open(s[1]) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip empty lines and comments
                    
                cols = line.split('\t')
                if len(cols) >= 9:  # Only include lines with at least 9 columns
                    split.append(cols)
                else:
                    if opts.v:
                        print("Warning: Skipping malformed line %d in %s (fewer than 9 columns): %s" % 
                              (line_num, s[1], line))

        ## i = numLine; v = corresponding i-th GTF row
        for i,v in enumerate(split):
            if is_transcript(v):
                try:
                    # Use the new getTranscriptID function
                    t_id = getTranscriptID(v[8])
                    if not t_id:
                        raise Exception("Could not find transcript_id")
                    g_id = getGeneID(v[8], v[0], t_id)
                except Exception as e:
                    print("Problem parsing file %s at line %d:\n%s\n" % (s[1], i+1, v))
                    print("Error: ", str(e))
                    continue  # Skip this problematic transcript instead of exiting
                    
                geneIDs[t_id]=g_id
                if not RE_STRING.match(g_id):
                    try:
                        badGenes.append([v[0],v[6], t_id, g_id, min(int(v[3]),int(v[4])), max(int(v[3]),int(v[4]))]) #chromosome, strand, cluster/transcript id, start, end
                        j=i+1
                        while j<len(split):
                            if len(split[j]) < 3:  # Check if this line has enough columns
                                j+=1
                                continue
                                
                            if split[j][2]=="exon":
                                try:
                                    badGenes[len(badGenes)-1].append((min(int(split[j][3]), int(split[j][4])), max(int(split[j][3]), int(split[j][4]))))
                                except (ValueError, IndexError) as e:
                                    if opts.v:
                                        print("Warning: Could not parse exon coordinates at line %d: %s" % (j+1, split[j]))
                                j+=1
                            else:
                                break  # Not an exon line, so break
                    except (ValueError, IndexError) as e:
                        if opts.v:
                            print("Warning: Could not process transcript at line %d: %s" % (i+1, v))
                            print("Error: ", str(e))

    except Exception as e:
        warnings.warn("Error processing file %s: %s" % (s[1], str(e)))
        continue  # Try next sample instead of exiting

    if len(badGenes) > 0:
        break  # We found the "bad" genes!

##THE CLUSTERING BEGINS!##
if opts.cluster and len(badGenes)>0:
    clusters=[] #lists of lists (could be sets) or something of transcripts
    badGenes.sort(key=itemgetter(3)) #sort by start coord...?
    i=0
    while i<len(badGenes): #rather un-pythonic
        temp_cluster=[badGenes[i]]

        k=0
        while k<len(temp_cluster):
            j=i+1
            while j<len(badGenes):
                if t_overlap(temp_cluster[k], badGenes[j]):
                    temp_cluster.append(badGenes[j])
                    del badGenes[j]
                else:
                    j+=1
            k+=1
        if len(temp_cluster)>1:
            clusters.append([t[2] for t in temp_cluster])
        i+=1

    if opts.v:
        print("Found %d clusters" % len(clusters))

    for c in clusters:
        c.sort()

    clusters.sort(key=itemgetter(0))
    legend=[]
    for u,c in enumerate(clusters):
        my_ID=opts.key+str((u+1))
        legend.append(list(itertools.chain.from_iterable([[my_ID],c]))) #my_ID, clustered transcript IDs
        for t in c:
            geneIDs[t]=my_ID
##            geneIDs[t]="|".join(c) #duct-tape transcript IDs together, disregarding ref_gene_names and things like that

    with open(opts.legend, 'w') as l_file:
        my_writer=csv.writer(l_file)
        my_writer.writerows(legend)

geneDict=defaultdict(lambda: defaultdict(lambda: 0)) #key=gene/cluster, value=dictionary with key=sample, value=summed counts
t_dict=defaultdict(lambda: defaultdict(lambda: 0))
guidesFile='' # file given with -G for the 1st sample

for q, s in enumerate(samples):
    if opts.v:
       print(">processing sample %s from file %s" % s)
    lno=0
    try:
        f = open(s[1])
        transcript_len=0
        t_id = None
        g_id = None
        coverage = None
        
        for l in f:
            lno+=1
            if l.startswith('#'):
                if lno==1:
                    ei=l.find('-e')
                    if ei<0:
                       print("Warning: sample file %s might not have been generated with -e option!" % (s[1]))
                       # Continue anyway instead of exiting
                    gf=RE_GFILE.search(l)
                    if gf:
                       gfile=gf.group(1)
                       if guidesFile:
                          if gfile != guidesFile:
                             print("Warning: sample file %s generated with a different -G file (%s) than the first sample (%s)" % (s[1], gfile, guidesFile))
                       else:
                          guidesFile=gfile
                    else:
                       print("Warning: sample %s might not have been processed with -G option!" % (s[1]))
                       # Continue anyway instead of exiting
                continue
            
            # Skip empty lines
            if not l.strip():
                continue
                
            v=l.split('\t')
            
            # Skip malformed lines
            if len(v) < 9:
                if opts.v:
                    print("Warning: Skipping malformed line %d in %s (fewer than 9 columns)" % (lno, s[1]))
                continue
                
            if v[2]=="transcript":
                if transcript_len>0 and t_id is not None and coverage is not None:
                    # Only process if we have valid transcript information
                    t_dict[t_id][s[0]] = int(ceil(coverage*transcript_len/read_len))
                
                # Reset values for the new transcript
                transcript_len = 0
                
                try:
                    # Use the getTranscriptID function for robust parsing
                    t_id = getTranscriptID(v[8])
                    if not t_id:
                        # One more attempt with a direct search
                        t_id = RE_TRANSCRIPT_ID.search(v[8]).group(1)
                except:
                    if opts.v:
                        print("Warning: Could not extract transcript_id from line %d: %s" % (lno, v[8]))
                    t_id = None
                    continue  # Skip this transcript
                
                try:
                    g_id = getGeneID(v[8], v[0], t_id)
                except:
                    if opts.v:
                        print("Warning: Could not extract gene_id from line %d: %s" % (lno, v[8]))
                    g_id = t_id  # Use transcript ID as fallback
                
                try:
                    coverage = getCov(v[8])
                except:
                    if opts.v:
                        print("Warning: Could not extract coverage from line %d: %s" % (lno, v[8]))
                    coverage = 0.0
                
            elif v[2]=="exon":
                # Only process exons if we have a valid transcript ID
                if t_id is not None:
                    try:
                        transcript_len += int(v[4])-int(v[3])+1  # because end coordinates are inclusive in GTF
                    except (ValueError, IndexError):
                        if opts.v:
                            print("Warning: Could not parse exon coordinates at line %d" % lno)

        # Process the last transcript in the file
        if transcript_len>0 and t_id is not None and coverage is not None:
            t_dict[t_id][s[0]] = int(ceil(coverage*transcript_len/read_len))

        f.close()

    except Exception as e:
        warnings.warn("Error processing file %s: %s" % (s[1], str(e)))
        continue  # Try next sample instead of exiting

    # Build the gene-level counts from transcript-level counts
    for i,v in t_dict.items():
        try:
            if i in geneIDs and s[0] in v:
                geneDict[geneIDs[i]][s[0]] += v[s[0]]
        except KeyError as e:
            if opts.v:
                print("Warning: Problem processing transcript %s: %s" % (i, str(e)))

if opts.v:
   print("..writing %s " % (opts.t))
   
try:
    with open(opts.t, 'w') as csvfile:
       my_writer = csv.DictWriter(csvfile, fieldnames = ["transcript_id"] + [x for x,y in samples])
       my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
       for i in t_dict:
            t_dict[i]["transcript_id"] = i
            my_writer.writerow(t_dict[i])
except Exception as e:
    print("Error writing transcript count matrix: %s" % str(e))
    sys.exit(1)
    
if opts.v:
   print("..writing %s " % (opts.g))
   
try:
    with open(opts.g, 'w') as csvfile:
       my_writer = csv.DictWriter(csvfile, fieldnames = ["gene_id"] + [x for x,y in samples])
       my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
       for i in geneDict:
            geneDict[i]["gene_id"] = i  # add gene_id to row
            my_writer.writerow(geneDict[i])
except Exception as e:
    print("Error writing gene count matrix: %s" % str(e))
    sys.exit(1)
    
if opts.v:
   print("All done.")