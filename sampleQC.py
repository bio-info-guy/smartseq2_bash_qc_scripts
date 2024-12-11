#!/usr/bin/python
import math
import sys
import random
from scipy.stats import variation
import scipy as sp
import argparse

def ArgParser():
    parser = argparse.ArgumentParser(description='SAM file and GTF parser for calculation of mapping qualities')
    parser.add_argument('-n', action="store", dest = "name", help='sample name')
    parser.add_argument('--sam', action="store", dest = "sam_path", help='SAM file path')
    parser.add_argument('--gtf', action="store", dest="gtf_path", help='GTF file path, must be same GTF used in mapping')
    parser.add_argument('--se', action="store_true", dest="single", default = False, help='option to indicate mapped reads are single ended')
    parser.add_argument('--rsem', action="store_true", dest="rsem", default=False, help='whether to use RSEM calculated expression for analysis' )
    parser.add_argument('--rsem_ct', action="store", dest="rsem_ct_path", help='RSEM read count file path Genes x Samples, must include sample name')
    parser.add_argument('--rsem_len', action="store", dest="rsem_len_path", help='RSEM gene length file path Genes x Samples, must include sample name')
    parser.add_argument('--rsem_tpm', action="store", dest="rsem_tpm_path", help='RSEM TPM file path Genes x Samples, must include sample name')
    parser.add_argument('--spike', action="store_true", dest="spike", default=False, help = 'this option tells the program to only look at spike in reads, only supports ERCC type spike-ins that are annotated in the GTF file')
    parser.add_argument('--sampling', action="store", dest="sampling", default=0, type=int, help = 'an integer indicating the amount of unique reads to sample to calculate quality, a value of 0 uses all reads')
    parser.add_argument('--gap', action="store", dest="gap_size", default=5, type=int, help = 'number of consecutive gene nt without coverage to be considered a gap')
    parser.add_argument('--window', action="store", dest="window", default=1, type=int, help = 'sliding window size to calculate evenness of coverage')
    results = parser.parse_args()
    if results.rsem:
        if results.rsem_ct_path is None or results.rsem_tpm_path is None or results.rsem_len_path is None:
            parser.print_help(sys.stderr)
            print('INVALID INPUT: please provide RSEM files if the rsem option is selected')
            sys.exit(0)
    return results
# cv #
# input: list of numbers (sum greater than zero)
# output: coeeficient of variation of the list of numbers

def cv(lis):
    mean = sum(lis)/float(len(lis))
    var = []
    variance = 0
    for num in lis:
        var.append((num-mean)**2)
    variance = sum(var)/(len(var))
    cffv = math.sqrt(variance)/float(mean)
    return cffv

# merge #
# merge function merges two tuples of (start, end) format if overlap exists
def merge(item1, item2):
    a = item1
    b = item2
    if a[0] <= b[0] and a[1] >= b[1]:
        return a
    elif a[0] < b[0] and a[1] < b[1]:
        if a[1] >= b[0]:
            new = (a[0], b[1])
            return new
        else:
            return 0
    elif a[0] >= b[0] and a[1] <= b[1]:
        return b
    elif a[0] > b[0] and a[1] > b[1]:
        if a[0] <= b[1]:
            new = (b[0],a[1])
            return new
        else:
            return 0


# exon #
# exon function merges all tuples of (start, end) format within a list to give an agregate list of (start, end) tuples
def exon(lis, use=False):
    cords = lis[:]
    f = True
    while f:
        f = False
        i = 0
        while i <= len(cords) - 2:
            j = i+1
            while j <= len(cords) - 1:
                new = merge(cords[i],cords[j])
                if new != 0:
                    cords[i]=new
                    cords.pop(j)
                    f = True
                    continue
                else:
                    j = j+1
            i = i+1
    m = []
    for item in cords:
        m.append(item[1]-item[0]+1)
    if use:
        return cords
    else:
        return m

   
# RSEM #
# read gene expression or gene length input from RSEM output
def RSEM(file, sample_name):
    first_line = open(file).readline().rstrip().split("\t")
    col = first_line.index(sample_name)
    attr = {}
    for line in open(file):
        current_gene = line.rstrip().split("\t")
        if sample_name in current_gene:
            continue
        else:
            attr[current_gene[0]] = float(current_gene[col])
    return attr


# GTF #
# takes a gtf file and stores its properties in python for usage
class GTF:
    def getexon(self, file):
        chromo = {}
        _genelen_ = {}
        _exon_ = {}#aggregate exons with overlapped regions merged
        _exon = {}#non-merged exons
        _pos_ = {}
        for line in open(file):
            if "#!" in line:
                continue
            elif line.split("\t")[2] != "exon":
                continue
            else:
                chrom = line.split("\t")[0]
                current_exon = line.split("\t")
                if chrom in chromo.keys():
                    gene_name = current_exon[8].split(";")[0].split('"')[1]
                    try:
                        chromo[chrom][gene_name].append((int(current_exon[3]), int(current_exon[4])))
                    except:
                        chromo[chrom][gene_name] = [(int(current_exon[3]), int(current_exon[4]))]
                else:
                    chromo[chrom]={}
                    gene_name = current_exon[8].split(";")[0].split('"')[1]
                    chromo[chrom][gene_name] = [(int(current_exon[3]), int(current_exon[4]))]

        for chro in chromo.keys():
            _pos_[chro] = {}
            for key in chromo[chro].keys():
                _exon[key] = (chro, chromo[chro][key])
                _exon_[key] = (chro, exon(chromo[chro][key], True))
                _genelen_[key] = exon(chromo[chro][key])[0]
        for key in _exon_.keys():
            for (start, end) in _exon_[key][1]:
                while start <= end:
                    try:
                        _pos_[_exon_[key][0]][start].append(key)
                    except:
                        _pos_[_exon_[key][0]][start] = [key]
                    start = start+1
        return (_exon_, _genelen_, _pos_, _exon)

    def __init__(self, value):
        self.data = self.getexon(value)

    def length(self, gene=False):
        if gene:
            return self.data[1][gene]
        else:
            return self.data[1]

    def exon(self, gene=False):
        if gene:
            return self.data[0][gene]
        else:
            return self.data[0]

    def position(self, pos=False):
        if pos:
            return self.data[2][pos]
        else:
            overlap = set() # returns an set of overlapped genes in the annotation file
            for chrom in self.data[2].keys():
                for item in self.data[2][chrom]:
                    if len(self.data[2][chrom][item]) >= 2 :
                        overlap.add(str(self.data[2][chrom][item]))
            return overlap


# SAM Class #
# takes a sam file and stores its properties in python for usage    
class SAM:
    # CIGAR function: for each line of a SAM file, take the start position of the read and the CIGAR string and converts into list of tuples denoting (start, end) segments.
    # example: start=1, cigar=10M15N100M would return [(1,10), (26,125)]
    def CIGAR(self, start, cigar):
        cigar_lis = list(cigar)
        ref = []
        ss = int(start)
        num = ""
        for char in cigar_lis:
            if char.isdigit():
                num = num+char
            else:
                if char in "SM":
                    ref.append((ss, ss+int(num)-1))
                    ss = ss+int(num)
                elif char in "ND":
                    ss = ss+int(num)
                num = ""
        return ref

    def sample(self, nam_set, num, bootstrap=True):
        read_nam = []
        size = num
        for nam in nam_set:
            read_nam.append(nam)
        if len(read_nam) < num and bootstrap:
            rand_records = sorted(random.choices(read_nam, k = size))
        elif len(read_nam) < num:
            size = len(nam_set)
            rand_records = sorted(random.sample(read_nam, k = size))
        else:
            rand_records = sorted(random.sample(read_nam, k = size))
        return rand_records, size

    def filter(self, line, single=False, spike_only = False):
        current_read = line.split("\t")
        if (spike_only and "ERCC" in line) or (not spike_only and "ERCC" not in line):
            consider = True
        else:
            consider = False
        if "YT:Z:U" in line and current_read[5] == "*" and not single:
            return False
        elif "ERCC" in line and not spike_only:
            return "ERCC"
        elif (current_read[6] == "=" or single) and consider:
            if not single:
                a = max(int(current_read[3]), int(current_read[7]))
                b = min(int(current_read[3]), int(current_read[7]))
                marker = current_read[2]+":"+str(a) + ":" + str(b)
            else:
                marker = current_read[2] + ":" + current_read[3]
            if "NH:i:1" in current_read and current_read[8] != "0":
                return (marker, "uniq")
            elif "NH:i:1" in current_read and int(current_read[4]) >= 60 and single:
                return (marker, "uniq")
            else:
                return marker
        elif consider:
            a = max(int(current_read[3]), int(current_read[7]))
            b = min(int(current_read[3]), int(current_read[7]))
            chromo = sorted([current_read[6], current_read[2]])
            if current_read[6] == "*":
                return current_read[2]+":"+current_read[3]
            else:
                return chromo[0] + ":" + chromo[1] + ":" + str(a) + ":" + str(b)
        else:
            return False
        
    def readSAM(self, file, single=False, spike_only=False):
        _allreads_ = {}
        _uniq_ = {}
        read = set()
        spike = set()
        nopass = set()
        for line in open(file):
            if line.startswith("@"):
                continue
            else:
                current_read = line.split("\t")
                read.add(current_read[0])
                type = self.filter(line.rstrip(), single, spike_only)
                if not type:
                    if spike_only:
                        continue
                    nopass.add(current_read[0])
                    _allreads_[current_read[0]] = set()
                elif "ERCC" in type and not spike_only:
                    spike.add(current_read[0])
                    _allreads_[current_read[0]] = set()
                else:
                    if "uniq" in type:
                        try:
                            _allreads_[current_read[0]].add(type[0])
                        except:
                            _allreads_[current_read[0]] = set()
                            _allreads_[current_read[0]].add(type[0])
                        try:
                            segments = self.CIGAR(current_read[3], current_read[5])
                            for seg in segments:
                                _uniq_[current_read[0]].append(seg)                                
                        except:
                            _uniq_[current_read[0]]=[current_read[2]]
                            segments = self.CIGAR(current_read[3], current_read[5])
                            for seg in segments:
                                _uniq_[current_read[0]].append(seg)
                    else:
                        try:
                            _allreads_[current_read[0]].add(type)
                        except:
                            _allreads_[current_read[0]] = set()
                            _allreads_[current_read[0]].add(type)
        for key in _uniq_.keys():
            if len(_uniq_[key]) > 2:
                _uniq_[key] = [_uniq_[key][0], exon(_uniq_[key][1:], True)]
            else:
                _uniq_[key] = [_uniq_[key][0], [_uniq_[key][1]]]
        if spike_only:
            read = set(_allreads_.keys())
        return (read, _allreads_, _uniq_, spike, nopass)

    def __init__(self, value, single=False, spike=False):
        self.data = self.readSAM(value, single, spike)


    def num_read(self):
        return len(self.data[0])

    def spike_rate(self):
        return round(len(self.data[3])/float(len(self.data[0])), 5)

    def uniq_rate(self):
        return round(len(self.data[2].keys())/float(len(self.data[0])-len(self.data[3])), 4)

    def display(self):
        print(self.data[0])
        print(self.data[1])
        for item in self.data[2].keys():
            if len(self.data[2][item][1]) > 1:
                print(item)
                print(self.data[2][item])

    def complexity(self, num, single=False):
        (smpl_read, size)=self.sample(self.data[0], num, bootstrap=False)
        uniq_read_pos = set()
#       all_uniq=set()
        j = 0
        pos = set()
        for key in self.data[1].keys():
            if len(self.data[1][key]) == 1:
                j = j+1
        for key in smpl_read:
            try:
                for item in self.data[1][key]:
                    pos.add(item)
                    if len(self.data[1][key]) == 1:
                        if len(item.split(":")) == 3:
                            uniq_read_pos.add(item)
                        elif len(item.split(":")) == 2 and single:
                            uniq_read_pos.add(item)
#                        all_uniq.add(item)
            except:
                continue
        complex = round(len(uniq_read_pos)/float(size), 3)
        return complex


# GeneCount #
# input: GTF class and SAM class instance
# output: dictionary of genes and their respective read counts
# NOTE:
# union: htseq-count union method
# intersect: htseq-count intersect method
# non-intersect:read count of read that map to the overlap of two genes are proportionally divided between the two genes based the ratio of two gene's nonoverlapped counts
# slight modification can make the output same as HTSeq-count output by turning off the overlap counting mechanism
def GeneCount(GTF, SAM, ct_opt="union"):
    reads = SAM.data[2]
    gtf = GTF.data[2]
    _count_ = {}
    _overlap_ = {}
    for gene in GTF.data[1].keys():
        _count_[gene] = 0
    for key in reads.keys():
        # step 0: for each read in uniq_reads, retreive the min position value and set as st, the max position value and set as en 
        chrom = reads[key][0]
        segs = []
        for (start, end) in reads[key][1]:
            segs.append(start)
            segs.append(end)
        st = min(segs)
        en = max(segs)
        try:
            st_gene = gtf[chrom][st]
            en_gene = gtf[chrom][en]
            #case 1:
            if len(st_gene) == 1 and len(en_gene) == 1:
                if st_gene[0] == en_gene[0]:
                    try:
                        _count_[st_gene[0]] += 1
                    except:
                        _count_[st_gene[0]] = 1
                else:
                    try:
                        _overlap_[st_gene[0] + ":" + en_gene[0]] += 1
                    except:
                        _overlap_[st_gene[0] + ":" + en_gene[0]] = 1
            elif len(st_gene) == 1 and len(en_gene) > 1 and ct_opt == "intersect":
                try:
                    _count_[st_gene[0]] += 1
                except:
                    _count_[st_gene[0]] = 1
            elif len(en_gene) == 1 and len(st_gene) > 1 and ct_opt == "intersect":
                try:
                    _count_[en_gene[0]] += 1
                except:
                    _count_[en_gene[0]] = 1
            elif ct_opt == "non-intersect":
                try:
                    _overlap_[st_gene[0] + ":" + st_gene[1]] += 1
                except:
                    _overlap_[st_gene[0] + ":" + st_gene[1]] = 1
        except:
            try:
                st_gene = gtf[chrom][st]
                if len(st_gene) == 1:
                    try:
                        _count_[st_gene[0]] += 1
                    except:
                        _count_[st_gene[0]] = 1
                elif ct_opt == "non-intersect":
                    try:
                        _overlap_[st_gene[0]+":"+st_gene[1]] += 1
                    except:
                        _overlap_[st_gene[0]+":"+st_gene[1]] = 1
            except:
                try:
                    en_gene = gtf[chrom][en]
                    if len(en_gene) == 1:
                        try:
                            _count_[en_gene[0]] += 1
                        except:
                            _count_[en_gene[0]] = 1
                    elif ct_opt == "non-intersect":
                        try:
                            _overlap_[en_gene[0]+":"+en_gene[1]] += 1
                        except:
                            _overlap_[en_gene[0]+":"+en_gene[1]] = 1
                except:
                    continue
    if ct_opt == "non-intersect":
        for key in _overlap_.keys():
            gene1 = key.split(":")[0]
            gene2 = key.split(":")[1]
            try:
                ratio = _count_[gene1]/float(_count_[gene2] + _count_[gene1])
                _count_[gene1] = int(_count_[gene1] + math.floor(_overlap_[key] * ratio))
                _count_[gene2] = int(_count_[gene2] + math.floor(_overlap_[key] * (1-ratio)))
            except:
                continue
    return _count_


# poswise depth function #
# input SAM class object, output dictionary pos_dp: keys are chromosomes, each corresponding item of key is dictionary of position depth
# thus pos_dp[chrI][1000] would be the coverage depth for the sample at chrI 1000 on the genome.
def poswise_depth(SAM, sampling = 0):
    uniq = SAM.data[2]
    if sampling:
        assert sampling >= 1000, 'please sample more reads than 1000, else sampled reads are not reliable'
        reads, _ = SAM.sample(SAM.data[0], sampling)
        
    else:
        reads = uniq.keys()
    pos_dp = {}
    for key in reads:
        try:
            chromo = uniq[key][0]
            segs = uniq[key][1]
        except:
            continue
        try:
            for st, en in segs:
                while st <= en:
                    try:
                        pos_dp[chromo][st] += 1
                    except:
                        pos_dp[chromo][st] = 1
                    st += 1
        except:
            pos_dp[chromo] = {}
            for st, en in segs:
                while st <= en:
                    try:
                        pos_dp[chromo][st] += 1
                    except:
                        pos_dp[chromo][st] = 1
                    st += 1
    return pos_dp


# rpkm function #
# input GTF class object, and gene count dictionary obtained from GeneCount function
# output dictionary norm_count: keys are genes, each corresponding item of key is the rpkm normalized read values
# norm_count[gene1] would thus be the normalized read count of gene1 for sample
def normalize(count, genlen, tpm=True):
    if not tpm:
        total = sum(count[key] for key in count.keys())
    else:
        total = sum([count[key]/float(genlen[key]) for key in count.keys() if genlen[key] != 0])
    norm_count = {}
    for key in count.keys():
        if genlen[key] == 0:
            norm_count[key] = 0
        else:
            if not tpm:
                norm_count[key] = round(count[key]*10**9/(float(genlen[key]*total)), 4)
            else:
                norm_count[key] = round((count[key]/genlen[key])*10**6/total, 4)
    return norm_count


# coverage function #
# usage: return coverage information of a gene
# input consists of:
    # chromo -> chromosome gene is located
    # exon_list -> aggregate list of the gene exon's
    # cov_dict -> coverage dictionary returned by the poswise_depth function
    # window -> window of coverage, a number of 10 would mean look at every tenth base starting from the first in exon_list
# output: coverage is a list containing the coverage information of the gene
def coverage(chromo, exon_list, cov_dict, window):
    # return coverage information for a certain gene given a fixed window of bases to look at
    num = int(window)
    coverage = []
    for (head, tail) in exon_list:
        cov = []
        pos = head
        while pos <= tail:
            try:
                cov.append(cov_dict[chromo][pos])
            except:
                cov.append(0)
            pos = pos+num
        if sum(cov) == 0:
            continue
        else:
            coverage.append(variation(sp.array(cov)))
    if len(coverage) == 0:
        return 0
    else:
        return sum(coverage)/float(len(coverage))


# gaps #
# calculate the number of gaps within a gene's exons
# usage: return number of uncovered continuous positions within a gene that are longer than length=gap
# input consists of:                                                                                                                                                              
    # chromo -> chromosome gene is located                                                                                                                                         
    # exon_list -> aggregate list of the gene exon's                                                                                                                               
    # cov_dict -> coverage dictionary returned by the poswise_depth function                                                                                                       
    # gap -> gap length, a gap length of 5 would mean if a continuous  region of longer than 5bp in the gene was not covered, it would be counted as a gap
# output: gap_num is  number of gaps of the gene
def gaps(chromo, exon_list, cov_dict, gap):
    num = int(gap)
    gap_num = 0
    j = 0
    translength = 0
    for (head, tail) in exon_list:
        counter = 0
        pos = head
        while pos <= tail:
            if counter == 5:
                translength = translength+(tail-head+1) 
            new = False
            try:
                if cov_dict[chromo][pos] > 0:
                    counter +=1
                    j = 0
                    pos += 1
                    new = False
                    continue
            except:
                j = j+1
                if j == num:
                    new = True
                else:
                    new = False
            if j >= num and new == True:
                gap_num = gap_num + 1
            pos = pos+1
    return gap_num, translength


# even_gaps #
# Function for calculation of evenness of coverage and gap statistics of the BAM alignment file produced by HISAT2.
# The specific calculations of the statistics are shown below:
# COVERAGE STATISTIC:
    # the coverage statistics looks at for each gene the evenness of coverage at each position
    # By moving moving across a window of user defined length, it takes all the position depth and calculate the coefficient of variation for the gene:
    # i.e: if window=10, and gene body spans chr1 0-100, it will look along the gene body at 0, 10, 20, 30, etc, and record the coverage depth at each of the positions.
    # Finally, it will output the coefficient of variaition for these 10 numbers
# GAP STATISTIC:
    # the gap statistics looks at for each gene the number of gaps within its gene body, normalized by the gene body effective length.
    # a gap is counted each time a region of longer than a user predefined length contains 0 coverage, by default it is 5.
    # the gene body effective length can be either two choices:
    # if aggregate is TRUE:
    # the effective length is the total length of aggregate exons that contain coverage of at least 5 positions.
    # else:
    # the effective length is the input from the output of software such as RSEM or SALFISH (RECOMMENDED)
def even_gaps(GTF, SAM, gene_count, gap, window, rsem=True, rsem_tpm=None, rsem_len=None, sampling = 0, spike=False):
    len_gap = int(gap)
    merged_exons = GTF.data[0]
    exons = GTF.data[3]
    pos_dict = poswise_depth(SAM, sampling = sampling)
    # calculate gaps for each gene, and store in dictionary
    gene_gaps = {}
    gene_cov = {}
    exonlen = {}
    if not rsem:
        out = open(SAMPLE_NAME + ".len", "w")
        for key in exons.keys():
            gene_gaps[key], exonlen[key] = gaps(merged_exons[key][0], merged_exons[key][1], pos_dict, len_gap) # retrieve the gene effective transcription length 
            gene_cov[key] = coverage(exons[key][0], exons[key][1], pos_dict, window=1)
            out.write(key+"\t"+str(exonlen[key])+"\n")
        out.close()
        gene_expr = normalize(gene_count, exonlen) 
    else:
        for key in exons.keys():
            gene_gaps[key], _ = gaps(merged_exons[key][0], merged_exons[key][1], pos_dict, len_gap)
            gene_cov[key] = coverage(exons[key][0], exons[key][1], pos_dict, window=1)
        exonlen = rsem_len
        if rsem_tpm is None:
            assert rsem_len != None, 'please provide rsem_len as argument to even_gaps function'
            gene_expr = normalize(gene_count, exonlen)
        else:
            gene_expr = rsem_tpm
    if spike:
        for k in gene_expr:
            if "ERCC" not in k:
                gene_expr[k] = 0
    # put number of gaps into python list for averaging
    mean_gaps = round(sum([gene_gaps[gene]*gene_expr[gene] for gene in gene_gaps.keys()])/float(sum([gene_expr[gene] for gene in gene_gaps.keys()])), 4)
    # take the top half of most expressed genes and calculate average coefficient of variation of coverage
    top_genes = sorted(gene_expr, key = gene_expr.__getitem__, reverse=True)
    half = int(len([gene_expr[key] for key in gene_expr.keys() if gene_expr[key] != 0])/2) if not spike else len([gene_expr[key] for key in gene_expr.keys() if gene_expr[key] != 0])
    even = []
    for gene in top_genes[:half]:
        even.append(gene_cov[gene])
    even = round(sum(even)/half, 4)
    return even, mean_gaps


def main(sample_name, sam_file, gtf_file, single=False, rsem=True, ct_path=None, len_path=None, tpm_path=None, spike = False, sampling = 0, gap = 5, window = 1):
    sam = SAM(sam_file, single, spike = spike)
    gtf = GTF(gtf_file)
    if rsem:
        count = RSEM(ct_path, sample_name)
        len_file = RSEM(len_path, sample_name)
        tpm=RSEM(tpm_path, sample_name)
    else:
        count = GeneCount(gtf, sam, "intersect")
        len_file = None
        out_count = open(sample_name+".ct", "w")
        out_count.write("\t"+sample_name+"\n")
        for key in sorted(count.keys()):
            out_count.write(key+"\t"+str(count[key])+"\n")
        out_count.close()
    sensitivity = sum([1 for key in count.keys() if count[key] != 0])
    rounds = 10 if not spike else 50
    even = 0; gaps = 0
    for i in range(rounds):
        (e, g) = even_gaps(gtf, sam, count, gap, window, rsem=rsem, rsem_tpm = tpm, rsem_len=len_file, sampling=sampling, spike=spike)
        even += e*float(1/rounds); gaps += g*float(1/rounds);
    file_name = sample_name+'_sampling_'+str(sampling)+'_gap_'+str(gap)+'_window_'+str(window)
    if rsem:
        file_name=file_name+'_rsem'
    if spike:
        file_name=file_name+'_spike'
    if single:
        file_name = file_name + '_se'
    out_stat = open(file_name + ".stat", "w")
    complexity_sample = 1000000 if not spike else sampling
    if spike:
        out_stat.write(sample_name + "\t" +
                   str(sam.num_read()) + "\t" +
                   str(sam.uniq_rate()) + "\t" +
                   str(sensitivity)+"\t"+
                   str(sam.complexity(complexity_sample, single)) + "\t" +
                   str(round(even,4))+ "\t" +
                   str(round(gaps,4))+"\n")
    else:
         out_stat.write(sample_name + "\t" +
                   str(sam.num_read()) + "\t" +
                   str(sam.spike_rate()) + "\t" +
	           str(sam.uniq_rate()) + "\t" +
                   str(sensitivity)+"\t"+
                   str(sam.complexity(complexity_sample, single)) + "\t" +
                   str(round(even,4))+ "\t" +
                   str(round(gaps,4))+"\n")
    out_stat.close()


if __name__ == "__main__":
    results = ArgParser()
    SAMPLE_NAME = results.name
    main(sample_name = results.name,
         sam_file = results.sam_path,
         gtf_file = results.gtf_path,
         single = results.single,
         rsem = results.rsem,
         ct_path = results.rsem_ct_path,
         tpm_path = results.rsem_tpm_path,
         len_path = results.rsem_len_path,
         spike = results.spike,
         sampling = results.sampling,
         gap = results.gap_size,
         window = results.window)
