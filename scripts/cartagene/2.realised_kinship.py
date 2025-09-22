# Original script by Sharon R. Browning: https://faculty.washington.edu/sguy/ibd_relatedness.html
# Updated using Python 3's 2to3 tool

# Input ibd files should follow the current Beagle/RefinedIBD format which includes the ninth column with IBD length.
# The map file should follow the four column plink format (as for use in Beagle/RefinedIBD).
# The map file is used to obtain the full genome length, assuming each chromosome starts at the smallest starting position among the IBD segments, and ends at the largest ending position among the IBD segments.
# If included, the vcf file header is used to read in the IDs, so that the matrix is ordered according to that. Otherwise, the IDs present in the IBD file will be used, and the order will be somewhat arbitrary. Vcf header can be gzipped (and it's okay if it is the full vcf file).

# Outputs non-zero relatedess only if there is a segment of length at least cmthresh centiMorgans. (Use zero for this threshold to turn this feature off.)

# By default, outputs relatedness values (2*kinship), but can output kinship instead by including the keyword in the commandline.
# By default, outputs a square matrix of relatedness values (1's on diagonal for relatedness, 0.5 for kinship), but can output a lower triangular matrix instead by including the keyword in the commandline.

# Whether outputting square matrix or lower triangular, will output a header line and first column with IDs. By default, the header line starts with a filler entry "ID", to ensure that the ids in the header row line up with the corresponding entries. However, with the "nofill" option, this entry is omitted and the header line contains one fewer entry than the other lines (for the full matrix). This makes it easy to read the data into R with read.table. 

# usage: zcat myibdfile.ibd.gz | python relatedness_v1.py mapfile cmthresh [vcfhead=myfile.vcf[.gz]] [kinship] [lower] [nofill] > outfile

# The optional arguments can be included in any order after the required arguments.

import sys, gzip

# function for processing keywords from the commandline
def process_keyword(keyword,kinship,lower,vcffilename,fillhead):   
    if keyword=="kinship":
        kinship = True
    elif keyword=="lower":
        lower = True
    elif keyword[:8]=="vcfhead=":
        vcffilename = keyword[8:]
    elif keyword=="nofill":
        fillhead = False
    else:
        print("Error: Unrecognized keyword in commandline",keyword, file=sys.stderr)
        raise SystemExit
    return (kinship,lower,vcffilename,fillhead)

# process the command line
mapfile = open(sys.argv[1])
cmthresh = float(sys.argv[2])
kinship=False; lower=False; vcffilename=""; fillhead=True
for i in range(3,len(sys.argv)):
    (kinship,lower,vcffilename,fillhead)=process_keyword(sys.argv[i],kinship,lower,vcffilename,fillhead)

#print >> sys.stderr, kinship,lower,vcffilename,fillhead
    
# read in the IDs, if vcf header file is available
if vcffilename != "":
    if vcffilename[-3:]==".gz":
        vcffile = gzip.open(vcffilename)
    else:
        vcffile = open(vcffilename)
    for line in vcffile:
        if line[0]=="#" and line[1:6]=="CHROM":
            ids = line.split()[9:]
            break
    else:
        print("Error: Unable to find header line in file",vcffilename, file=sys.stderr)
        raise SystemExit
    
# read in the map file
gpos = dict()
for line in mapfile:
    bits = line.split()
    chrom = bits[0]
    if chrom not in gpos:
        gpos[chrom] = dict()
    gpos[chrom][int(bits[3])] = float(bits[2])

# function for converting bp to cM    
def interpolate(bp,chromgpos):
    set1 = [x for x in chromgpos if x < bp]
    set2 = [x for x in chromgpos if x > bp]
    if len(set1)>0 and len(set2)>0:
        bp1 = max(set1); bp2 = min(set2)
    elif len(set1)>0:
        bp2 = max(set1); set1.remove(bp2); bp1 = max(set1)
    else:
        bp1 = min(set2); set2.remove(bp1); bp2 = min(set2)
    # interpolate
    if bp1 == bp2: return chromgpos[bp1]
    else: return chromgpos[bp1]+(chromgpos[bp2]-chromgpos[bp1])*(bp - bp1)/float(bp2-bp1)

# start reading in the total IBD length (in cM), and also keep track of chromosome start and endpoints (in bp)
ibd = dict()
allids = set()
chromminbp = dict()
chrommaxbp = dict()
for line in sys.stdin:
        bits = line.split()
        chrom = bits[4]
        length = float(bits[8])
        id1 = bits[0]; id2 = bits[2]
        start = int(bits[5]); end = int(bits[6])
        if chrom not in chromminbp:
            chromminbp[chrom] = start
            chrommaxbp[chrom] = end
        else:
            chromminbp[chrom] = min(start,chromminbp[chrom])
            chrommaxbp[chrom] = max(end,chrommaxbp[chrom])
        if id1 not in ibd:
            ibd[id1] = dict()
            allids.add(id1)
        if id2 not in ibd[id1]:
            ibd[id1][id2] = []
            allids.add(id2)
        ibd[id1][id2].append(length)

if vcffilename != "":
    # check whether ids from header match those found in ibd file
    idsnotinibd = [x for x in ids if x not in allids]
    ibdnotinids = [x for x in allids if x not in ids]
    if len(idsnotinibd)>0:
        print(len(idsnotinibd), "IDs in vcf header are not in the IBD file. These IDs will necessarily have kinship of zero with all other individuals.", file=sys.stderr)
        if len(idsnotinibd)<5:
            for x in idsnotinibd:
                print(x, file=sys.stderr)
    if len(ibdnotinids)>0:
        print(len(ibdnotinids), "individuals in the IBD file are not included in the vcf header. These individuals will not be included in the output.", file=sys.stderr)
        if len(ibdnotinids)<5:
            for x in ibdnotinids:
                print(x, file=sys.stderr)
else:
    ids = list(allids)

totalchromlen = 0
for chrom in chromminbp:
    mincm = interpolate(chromminbp[chrom],gpos[chrom])
    maxcm = interpolate(chrommaxbp[chrom],gpos[chrom])
    totalchromlen += maxcm - mincm
if totalchromlen == 0:
    print("Error: Inferred genome length is zero (probably due to no IBD in input", file=sys.stderr)
    raise SystemExit

if kinship:
    denom = totalchromlen*4.0
else:
    denom = totalchromlen*2.0

# print header row
if fillhead:
    print('\t'.join(["ID"]+ids))
else:
    print('\t'.join(ids))
# print other rows
for i,id1 in enumerate(ids):
    myrow = [id1]
    if lower:
        jvals = list(range(i))
    else:
        jvals = list(range(len(ids)))
    for j in jvals:
        id2 = ids[j]
        thisibd = []
        if id1 in ibd and id2 in ibd[id1]: thisibd.extend(ibd[id1][id2])
        if id2 in ibd and id1 in ibd[id2]: thisibd.extend(ibd[id2][id1])
        # decide whether to include pair based on longest ibd segment
        if len(thisibd)==0 or max(thisibd) < cmthresh:
            num = 0
        else:
            num = sum(thisibd)
        if i==j:
            num = 2.0*totalchromlen
        myrow.append("%.3g"%(num/denom))
    print('\t'.join(myrow))
