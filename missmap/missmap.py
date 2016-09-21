# TODO:
# class close function
# post-alignment processing
# pair-end support?

import sys
import os
from subprocess import Popen, PIPE


# ============================================================================


# ============================================================================

def mesg(message, ofile=sys.stderr, fatal=False):
    print(message, flush=True, file=ofile)
    if fatal:
        sys.exit(0)

# ============================================================================

def exec_cmd(cmd, dry=False):
    mesg(':: ' + cmd)
    if dry:
        return
    ecode = os.system(cmd)
    if ecode:
        mesg('*FAILED* exit code={}'.format(ecode), fatal=True)

# ============================================================================

# define fasta class and variant class
# fasta class has ...
# for variant class the underlying repr is unique, but can have multiple interface

class Variant(object):
    """Class of one genetic variant.

    TODO:
        check parsimony,
        normalization,
        decomposition,
        unique representation

    Assuming variants on usually on forward strand.
    """

    def __init__(self, chrom, pos, ref, alt):
        # pos is 1-based int
        self.chrom = chrom
        self.pstart = pos-1
        self.pend = self.pstart + len(ref)
        self.ref = ref
        self.alt = alt
        if len(ref) == 1 and len(alt) == 1:
            self.vtype = 'SNP'
        elif len(ref) == 1 and len(alt) > 1 and ref == alt[0]:
            # if len(alt) < 1kbp, otherwise, SV
            self.vtype = 'INS'
        elif len(alt) == 1 and len(ref) > 1 and alt == ref[0]:
            # if len(ref) < 1kbp, otherwise, SV
            self.vtype = 'DEL'
        else:
            # MNP complex (e.g. SNP+INDEL)
            self.vtype = 'MIX'

    def __str__(self):
        # chrom:pos:ref:alt:type
        repr_str = ('{}:{}:{}:{}:{}'.format(
            self.chrom, self.pstart+1, self.ref, self.alt, self.vtype))
        return repr_str

# ============================================================================


class Fasta(object):

    def __init__(self, chrom=None, pos=0, seq=[], rev=False, compl=False):
        """Initialize a fasta object.

        pos is 0-based starting position
        seq is a list
        """
        self.chrom = chrom
        self.pos = pos
        self.seq = seq
        self.rev = rev
        self.compl = compl

    # -------------------------------------------------------------------------

    def __str__(self):

        repr_str = ('> {}:{}-{} rev={} compl={}\n{}'.format(
            self.chrom, self.pos+1, self.pos+len(self.seq),
            self.rev, self.compl, ''.join(self.seq)))
        return repr_str

    # -------------------------------------------------------------------------

    def load_ref(self, ref_fasta_name, chrom, start, end):
        """Get sequence in a region from fasta file.

        query region: [start, end) 0-based half-open interval.
        """

        if start >= end:
            mesg('*ERROR* Variant.load_ref() wrong query range {}-{} '
                 '(start>=end)'.format(start, end), fatal=True)

        # open fai and get offset of region
        chrom_found = False
        for line in open(ref_fasta_name + '.fai'):
            tchrom, tlen, toffset, tnbase, tnchar = line.rstrip().split('\t')
            tlen, toffset, tnbase, tnchar = (int(tlen), int(toffset),
                                             int(tnbase), int(tnchar))
            if tchrom == chrom:
                chrom_found = True
                break

        # check if chrom in fasta
        if not chrom_found:
            mesg('*ERROR* Variant.load_ref(): chrom {} is not found in '
                 'fasta {}'.format(chrom, ref_fasta_name), fatal=True)

        # check if range out of chromosome boundary
        if start < 0:
            start = 0
            mesg('*WARNING* Fasta.load_ref(): start position is less than 0.')
        if end > tlen - 1:
            end = tlen - 1
            mesg('*WARNING* Fasta.load_ref(): end position is beyond chrom '
                 'length.')

        # read fasta and fetch sequence
        # calculate byte-offset of region in the reference fasta
        tnextra = tnchar - tnbase
        offset_start = toffset + start + tnextra*(start//tnbase)
        offset_end = toffset + (end-1) + tnextra*((end-1)//tnbase)

        seq_str = ''
        with open(ref_fasta_name) as f:
            f.seek(offset_start)
            seq_str = f.read(offset_end-offset_start+1).replace('\n','')

        # assign results
        self.chrom = chrom
        self.pos = start
        self.seq = list(seq_str)

    # -------------------------------------------------------------------------

    def rev_compl(self, rev=True, compl=True):
        """Return the reverse/complement sequence."""
        pair = {
            'A':'T', 'T':'A', 'C':'G', 'G':'C',
            'a':'t', 't':'a', 'c':'g', 'g':'c',
        }
        seq_rev = []

        tseq = self.seq[::-1] if rev else self.seq
        for s in tseq:
            if rev:
                s = s[::-1]
            ss = ''
            for c in s:
                ss += pair[c] if c in pair and compl else c
            seq_rev.append(ss)
        new_fasta = Fasta(self.chrom, self.pos, seq_rev, rev=rev, compl=compl)
        return new_fasta

    # -------------------------------------------------------------------------

    def add_var(self, var):
        """Add a variant to fasta sequence.

        var is a Variant object
        change in Fasta.seq:
            replace reference base(s) with ""
            replace base at pos with alternate base(s) ***IN LOWER CASE***

        limitation:
            does not support add variant to reverse complement sequence
        """

        if not self.seq:
            mesg('*ERROR* fail to add variant to emtpy sequence', fatal=True)

        if self.rev or self.compl:
            mesg('*ERROR* adding variant to reverse complement strand is not '
                 'yet supported', fatal=True)

        if (var.pstart >= self.pos + len(self.seq) or
            var.pend < self.pos):
            mesg('*ERROR* Fasta.add_var() variant is out side of sequence',
                 fatal=True)

        for p in range(var.pstart, var.pend):
            self.seq[p-self.pos] = ""
        self.seq[var.pstart-self.pos] = var.alt.lower()


# ============================================================================
# ============================================================================
# ============================================================================

class Simulation(object):

    # -------------------------------------------------------------------------

    def __init__(self, read_len):
        # [pstart, pend) 0-based half-open
        self.var_list = []
        self.chrom = None
        self.pstart = 1e12
        self.pend = 0
        self.read_len = read_len
        self.fasta = None
        self.read_list = []

    # -------------------------------------------------------------------------

    def load_var_list(self, var_list_csv):
        with open(var_list_csv) as f:
            for line in f:
                data = line.rstrip().replace('"', '').split(',')
                chrom, pos, foo, ref, alt = data[:5]
                if not self.chrom:
                    self.chrom = chrom
                elif self.chrom != chrom:
                    mesg('*ERROR*: only support variants on the same chrom'
                         ' (for now).', fatal=True)

                var = Variant(chrom, int(pos), ref, alt)
                self.var_list.append(var)
                if var.pstart < self.pstart:
                    self.pstart = var.pstart
                if var.pend > self.pend:
                    self.pend = var.pend

    # -------------------------------------------------------------------------

    def simulate_reads(self, ref_fasta_name):

        self.fasta = Fasta()
        self.fasta.load_ref(ref_fasta_name, self.chrom, self.pstart, self.pend)

        for var in self.var_list:
            self.fasta.add_var(var)

        fwd_seq = ''.join(self.fasta.seq)
        # bwd_seq = ''.join(self.fasta.seq.rev_compl())

        # the position that it should align to??
        for i in range(0, len(forward_str)-self.read_len):
            read_list.append(fwd_seq[i:i+read_len])
            # read_list.append(bwd_seq[i:i+read_len])

    # -------------------------------------------------------------------------




# ============================================================================

class Alignment(object):

    # for the moment single end aligner

    def __init__(self, ref_fasta, aligner, exec_aligner,
                 nthread=4,
                 exec_samtools=None):
        """Initialize an alignment.

        aligner: string BWA-MEM/BWA-ALN
        """

        self.aligner = aligner
        self.config = {
            'exec_aligner': exec_aligner,
            'ref_fasta': ref_fasta,
            'nthread': nthread,
            'pair_end': False,
            'exec_samtools': 'samtools',  # assume in $PATH
            'post_proc': False,
        }
        if exec_samtools:
            self.config['exec_samtools'] = exec_samtools

    # -------------------------------------------------------------------------

    def make_rgstr(self, fq_in, sample_id, platform='ILLUMINA'):
        """ Derive RG string from fastq file.

        @RG\tID:${rgid}\tPU:${rgpu}\tSM:${rgsm}\tPL:${rgpl}

        In fastq:
            @HWI-7001446:480:C6BH4ANXX:3:1101:2304:1985 1:N:0:CCTCCT
        the first column should contain 7 fields:
            instrument id: HWI-7001446
            run id (flowcell_id): 480
            flowcell barcode: C6BH4ANXX
            lane id: 3
            tile id: 1101
            cluster x coordinate: 2304
            cluster y coordinate: 1985
        RG string will be:
            ID = instrument id : run id : ... : lane id
            PU = instrument id : run id : ... : flowcell barcode : lane id

        see https://en.wikipedia.org/wiki/FASTQ_format
        """

        rgdata = {}
        if fq_in.endswith('gz'):
            cmd = 'zcat ' + fq_in
        else:
            cmd = 'cat ' + fq_in
        proc = Popen(cmd, stdout=PIPE, shell=True, universal_newlines=True)
        for line in proc.stdout:
            seq_id = line.rstrip().split()[0].split(':')
            break
        proc.kill()

        if len(seq_id) != 7:
            mesg('*WARNING* read_in in {} does not contain 7 fields'
                 ' (RG string may be wrong)\n {}'.format(fq_in, seq_id))

        rgdata['instr_id'] = seq_id[0].lstrip('@')
        rgdata['run_id'] = seq_id[1]
        rgdata['flowcell_barcode'] = seq_id[2]
        rgdata['lane_id'] = seq_id[3]

        rgdata['ID'] = ('{}.{}.{}'.format(rgdata['instr_id'],
                                          rgdata['run_id'],
                                          rgdata['lane_id']))
        rgdata['PU'] = ('{}.{}.{}.{}'.format(rgdata['instr_id'],
                                             rgdata['run_id'],
                                             rgdata['flowcell_barcode'],
                                             rgdata['lane_id']))
        rgdata['SM'] = sample_id
        rgdata['PL'] = 'ILLUMINA'

        rgstr = ('@RG\tID:{}\tPU:{}\tSM:{}\tPL:{}'
                 .format(rgdata['ID'], rgdata['PU'], rgdata['SM'],
                         rgdata['PL']))
        return rgstr

    # -------------------------------------------------------------------------

    def run(self, sample_id, fq1_in, bam_out, fq2_in=None, rgstr=None):
        # read_id can not be longer than 256 chars

        # make rgstr
        if not rgstr:
            rgstr = self.make_rgstr(fq1_in, sample_id)

        # intermediate files
        sam_tmp = os.path.splitext(bam_out)[0] + '.sam'
        if self.aligner == 'BWA-ALN':
            sai_tmp = os.path.splitext(bam_out)[0] + '.sai'
        bam_tmp = bam_out + '.tmp'

        cmd_list = []
        # alignment
        if self.aligner == 'BWA-MEM':
            if not self.config['pair_end']:
                cmd = ('{} mem -t{} -M -R "{}" {} {} > {}'
                       .format(self.config['exec_aligner'],
                               self.config['nthread'], rgstr,
                               self.config['ref_fasta'], fq1_in, sam_tmp))
                cmd_list.append(cmd)
            else:
                # not implemented
                pass

        elif self.aligner == 'BWA-ALN':
            if not self.config['pair_end']:
                cmd = ('{} aln -t{} {} {} > {}'
                       .format(self.config['exec_aligner'],
                               self.config['nthread'],
                               self.config['ref_fasta'], fq1_in, sai_tmp))
                cmd_list.append(cmd)
                cmd = ('{} samse -r "{}" {} {} {} > {}'
                       .format(self.config['exec_aligner'], rgstr,
                               self.config['ref_fasta'], sai_tmp, fq1_in,
                               sam_tmp))
                cmd_list.append(cmd)
            else:
                # not implemented
                pass

        # covert sam to bam
        cmd = ('cat {} | {} view -bS - > {}'
               .format(sam_tmp, self.config['exec_samtools'], bam_tmp))
        cmd_list.append(cmd)

        # sort bam
        cmd = ('{} sort -T foo -O bam {} > {}'
               .format(self.config['exec_samtools'], bam_tmp, bam_out))
        cmd_list.append(cmd)

        # index bam
        cmd = ('{} index {}'.format(self.config['exec_samtools'], bam_out))
        cmd_list.append(cmd)

        # if clean up
        # cmd = ('rm {} {}'.format(sam_tmp, bam_tmp))
        # if self.aligner == 'BWA-ALN':
        #     cmd += ' {}'.format(sai_tmp)
        # cmd_list.append(cmd)

        # run commands
        for cmd in cmd_list:
            exec_cmd(cmd, dry=False)


    pass



# ============================================================================

class Report(object):

    pass


# ============================================================================
# ============================================================================



def main():

    REF_FASTA = ('/stornext/snfs4/1000-gen/DATA/REF/human/GRCh37/'
                 'human.g1k.v37.fa')
    REF_FASTA_BWT = ('/stornext/snfs4/1000-gen/DATA/REF/human/GRCh37/'
                     'bwa/human.g1k.v37.fa')

    # ===================== #
    # test Alignment module
    # ===================== #

    # BWA = '/stornext/snfs4/1000-gen/zhuoyih/software/bin/bwa'
    # fasta_name = '211487.fastq.gz'
    # sample_id = '211487'
    # bam_mem_name = '211487.bwa-mem.bam'
    # bam_aln_name = '211487.bwa-aln.bam'

    # ------------ test bwa mem ------------
    # aln = Alignment(ref_fasta=REF_FASTA_BWT, aligner='BWA-MEM',
    #                 exec_aligner=BWA)
    # aln.run('211487', '211487.fastq.gz', '211487.bwa-mem.bam')

    # ------------ test bwa mem ------------
    # aln = Alignment(ref_fasta=REF_FASTA_BWT, aligner='BWA-ALN',
    #                 exec_aligner=BWA)
    # aln.run('211487', '211487.fastq.gz', '211487.bwa-aln.bam')

    # sys.exit(0)



    # ref_fasta = pysam.FastaFile(ref_fasta_name)
    var_list = [
        # '22:46615915:GGTGTTTGCGGCT:G',
        '19:45406457:C:A',
        '19:45406459:G:C',
        # '19:45406461:TTCCACCTCCACC:T',
        '19:45406461:T:TTCCACC',
        '19:45406468:T:A',
        '19:45406469:C:A',
    ]

    fasta = Fasta()
    fasta.load_ref(REF_FASTA, '22', 46615714, 46615880)
    print(fasta, len(fasta.seq))
    fasta_rev = fasta.rev_compl(rev=True,compl=True)
    print(fasta_rev)

    # var = Variant('22', 46615778, 'A', 'C')
    # var = Variant('22', 46615778, 'A', 'AGAGT')
    var = Variant('22', 46615778, 'ATGACATAGAAGA', 'A')
    # var = Variant('22', 46615778, 'ATGACATAGAAGA', 'CGTCAGTGCGTACGTTAGTCGAT')
    fasta.add_var(var)
    print(fasta.seq)
    print(fasta, len(fasta.seq))






    pass


if __name__ == '__main__':
    main()













sys.exit(0)



















def reverse_seq(seq_list, flip=True):
    pair = {
        'A':'T', 'T':'A', 'C':'G', 'G':'C',
        'a':'t', 't':'a', 'c':'g', 'g':'c',
    }
    out_seq_list = []
    for s in seq_list:
        s_new = []
        for c in s:
            if c not in pair:
                s_new.append(c)
            else:
                s_new.append(pair[c])
        s_new = ''.join(s_new)
        if flip:
            out_seq_list.append(s_new[::-1])
        else:
            out_seq_list.append(s_new)

    if flip:
        out_seq_list = out_seq_list[::-1]

    return out_seq_list

# ============================================================================


# variant key like chrom:pstart(0-based inclusive):pend(0-based exclusive):substitution

# ============================================================================

def make_alt_sequence(var_list, ref_fasta, read_length):
    # var format: chrom:pos:ref:alt (ref/alt can be multi-bases)
    # better input var_list:
    # - sorted by position
    # - no overlap
    # - on the same chromosome
    # - pairwise distance between adjacent region is less than read_length
    # how to generate BQ string? (given at pos, model or random)
    # insertion size for pair end

    # check start and end position of var
    pmin = 1e16
    pmax = 0
    tchrom = ''
    for var in var_list:
        chrom, pos, ref, alt = var.split(':')
        if tchrom == '':
            tchrom = chrom
        else:
            if tchrom != chrom:
                mesg('*ERROR* variant need to be on the same chromosome.')
                sys.exit(0)
        pstart = int(pos) - 1  # 0-based inclusive
        pend = pstart + len(ref)  # 0-based exclusive
        if pstart < pmin:
            pmin = pstart
        if pend > pmax:
            pmax = pend

    # if variants are well separated, give warning
    if pmax - pmin > read_length:
        mesg('*WARNING* far separated variants may not be on the same read!')

    # reads have at least 1-base overlapping with variable region
    prange = [pmin - (read_length - 1), pmax + (read_length - 1)]

    # get reference seq
    ref_seq = list(ref_fasta.fetch(reference=chrom, start=prange[0],
                                   end=prange[1]).upper())
    alt_seq = ref_seq.copy()

    # apply variants sequentially (overlapping variants are not mutable)
    for var in var_list:
        chrom, pos, ref, alt = var.split(':')
        pstart = int(pos) - 1  # 0-based inclusive
        pend = pstart + len(ref)  # 0-based exclusive
        # remove ref allele
        for i in range(pstart-prange[0], pend-prange[0]):
            ref_seq[i] = ref_seq[i].lower()
            alt_seq[i] = '-'
        alt_seq[pstart-prange[0]] = alt.lower()

    # create the reverse strand of alternate sequence
    alt_seq_rev = reverse_seq(alt_seq)


    return [alt_seq, alt_seq_rev]

# ============================================================================

def write_fastq(fastq_name, seq_data, read_length, bq_mean=40, bq_stdev=0):

    # simulate BQ string
    bq_list = []
    for i in range(len(seq_data)):
        bq_list.append(chr(bq_mean+33))
    bq_list_rev = bq_list[::-1]

    # simulate read id
    read_id = "@SEQ_ID"

    # fragment the seq and write to disk
    # fout = fopen(fastq_name, 'w')
    # for i in range(0, len(alt_seq)-read_length):


    # fout.write()

# ============================================================================

    # print()
    # print(''.join(alt_seq))
    # print(''.join(alt_seq_rev))

    # # alt seq length or ref seq length? otherwise the output read length may
    # # not be the same
    # for i in range(0, len(alt_seq)-read_length):
    #     print(' '*i + ''.join(alt_seq)[i:i+read_length])
    # for i in range(0, len(alt_seq_rev)-read_length):
    #     print(' '*i + ''.join(alt_seq_rev)[i:i+read_length])






    # alt_seq = list(ref_seq)
    # for i in range(read_length-1, read_length-1+pmax-pmin):
    #     alt_seq[i] = '-'
    # alt_seq[read_length-1] = alt
    # print(ref_seq)
    # print(''.join(alt_seq))

    # return seq_str


# def simulate_fastq():
