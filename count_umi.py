import sys
import os
import argparse
import collections
import logging
import gzip
from functools import partial
import pandas as pd
import itertools
from dedup_umi import edit_distance

class Record:
    """A record representing a :term:`fastq` formatted record.

    Attributes
    ----------
    identifier : string
       Sequence identifier
    seq : string
       Sequence
    quals : string
       String representation of quality scores.
    format : string
       Quality score format. Can be one of ``sanger``,
       ``phred33``, ``phred64`` or ``solexa``.

    """
    def __init__(self, identifier, seq, quals, entry_format=None):
        self.identifier, self.seq, self.quals, entry_format = (
            identifier, seq, quals, entry_format)

    def __str__(self):
        return "@%s\n%s\n+\n%s" % (self.identifier, self.seq, self.quals)

def error(message):
    '''log error message, see the :mod:`logging` module'''
    logging.error(message)
    raise ValueError("Failed. Check the log file")

def openFile(filename, mode="r", create_dir=False):

    _, ext = os.path.splitext(filename)

    if ext.lower() in (".gz", ".z"):
        if sys.version_info.major >= 3:
            if mode == "r":
                return gzip.open(filename, 'rt', encoding="ascii")
            elif mode == "w":
                return gzip.open(filename, 'wt', compresslevel=6, encoding="ascii")
            else:
                raise NotImplementedError(
                    "mode '{}' not implemented".format(mode))
        else:
            return gzip.open(filename, mode, compresslevel=6)
    else:
        return open(filename, mode)

def addBarcodesToIdentifier(read, barcode, sgrna):
    '''extract the identifier from a read and append the UMI and
    cell barcode before the first space'''

    read_id = read.identifier.split(" ")

    if sgrna == "":
        read_id[0] = read_id[0] + "_" + barcode
    else:
        read_id[0] = read_id[0] + "_" + barcode + "_" + sgrna

    identifier = " ".join(read_id)

    return identifier

def getUserDefinedBarcodes(whitelist_csv, getErrorCorrection=False):
    '''get valid sgRNAs from sgrna-gene-mapping file formated as: id,sgrna,gene
    :param whitelist_csv:
    :param getErrorCorrection:
    :return:
    '''
    sgrna_whitelist = []
    gene_whitelist = []

    with openFile(whitelist_csv, "r") as inf:

        for line in inf:
            if line.startswith('#'):
                continue

            line = line.strip().split(",")
            whitelist_seq = line[1]
            whitelist_gene = line[2]
            sgrna_whitelist.append(whitelist_seq)
            gene_whitelist.append(whitelist_gene)

    return sgrna_whitelist, gene_whitelist

def get_substr_slices(umi_length, idx_size):
    '''
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function
    '''
    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    slices = []
    for s in sub_sizes:
        slices.append((offset, offset + s))
        offset += s
    return slices

def build_substr_idx(umis, umi_length, min_edit):
    '''
    Build a dictionary of nearest neighbours using substrings, can be used
    to reduce the number of pairwise comparisons.
    '''
    substr_idx = collections.defaultdict(
        lambda: collections.defaultdict(set))
    slices = get_substr_slices(umi_length, min_edit + 1)
    for idx in slices:
        for u in umis:
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(u)
    return substr_idx

def iter_nearest_neighbours(umis, substr_idx):
    '''
    Added by Matt 06/05/17
    use substring dict to get (approximately) all the nearest neighbours to
    each in a set of umis.
    '''
    for i, u in enumerate(umis, 1):
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            u_sub = u[slice(*idx)]
            neighbours = neighbours.union(substr_map[u_sub])
        neighbours.difference_update(umis[:i])
        for nbr in neighbours:
            yield u, nbr

class ExtractFilterAndUpdate:
    ''' A functor which extracts barcodes from a read(s), filters the
    read(s) and updates the read(s). Keeps track of events in
    read_counts Counter
    '''

    def _extract_5prime(self, sequence, read=1):
        if read == 1:
            return (sequence[:self.pattern_length],
                    sequence[self.pattern_length:])
        elif read == 2:
            return (sequence[:self.pattern_length2],
                    sequence[self.pattern_length2:])

    def _joiner_5prime(self, sequence, sample):
        return sample + sequence

    def _getBarcodesString(self, read1, read2=None):

        if self.pattern:
            bc1, sequence1 = self.extract(read1.seq)
            umi = "".join([bc1[x] for x in self.umi_bases])
            cell = "".join([bc1[x] for x in self.cell_bases])
            sample1 = "".join([bc1[x] for x in self.bc_bases])
            new_seq = self.joiner(sequence1, sample1)

        else:
            cell, umi, new_seq = ("",)*3

        if self.pattern2:
            bc2, sequence2 = self.extract(read2.seq, read=2)
            umi2 = "".join([bc2[x] for x in self.umi_bases2])
            cell2 = "".join([bc2[x] for x in self.cell_bases2])
            sample2 = "".join([bc2[x] for x in self.bc_bases2])
            new_seq2 = self.joiner(sequence2, sample2)

            cell += cell2
            umi += umi2
        else:
            new_seq2 = ""

        return cell, umi, new_seq, new_seq2


    def _getCellBarcodeString(self, read1, read2=None):

        if self.pattern:
            bc1, sequence1 = self.extract(read1.seq)
            cell = "".join([bc1[x] for x in self.cell_bases])
        else:
            cell = ""

        if self.pattern2:
            bc2, sequence2 = self.extract(read2.seq, read=2)
            cell2 = "".join([bc2[x] for x in self.cell_bases2])

            cell += cell2

        return cell

    def filterCellBarcode(self, cell):
        '''Filter out cell barcodes not in the whitelist'''

        if cell not in self.cell_whitelist:
            self.read_counts['Filtered cell barcode'] += 1
            return None,None

        ii = self.cell_whitelist.index(cell) ###
        gene = self.gene_whitelist[ii]   ###

        return cell, gene

    def __init__(self,
                 method="string",
                 pattern=None,
                 pattern2=None,
                 extract_cell=False,
                 filter_sgrna=True,
                 ):

        self.read_counts = collections.Counter()
        self.method = "string"
        self.pattern = pattern
        self.pattern2 = pattern2
        self.extract_cell = False
        self.filter_sgrna = True


        self.cell_whitelist = None  # These will be updated if required
        self.gene_whitelist = None #

        # If the pattern is a string we can identify the position of
        # the cell and umi bases at instantiation
        if method == "string":

            self.extract = self._extract_5prime
            self.joiner = self._joiner_5prime

            if pattern:
                if len(pattern.replace("N", "").replace("X", "").replace("C", "")) > 0:
                    raise ValueError("barcode pattern (%s) should only contain "
                                     "N/X/C characters" % pattern)
                self.pattern_length = len(pattern)
                self.umi_bases = [x for x in range(len(pattern)) if pattern[x] is "N"]
                self.bc_bases = [x for x in range(len(pattern)) if pattern[x] is "X"]
                self.cell_bases = [x for x in range(len(pattern)) if pattern[x] is "C"]

            if pattern2:
                if len(pattern2.replace("N", "").replace("X", "").replace("C", "")) > 0:
                    raise ValueError("barcode pattern2 (%s) should only contain "
                                     "N/X/C characters" % pattern2)
                self.pattern_length2 = len(pattern2)
                self.umi_bases2 = [x for x in range(len(pattern2)) if pattern2[x] is "N"]
                self.bc_bases2 = [x for x in range(len(pattern2)) if pattern2[x] is "X"]
                self.cell_bases2 = [x for x in range(len(pattern2)) if pattern2[x] is "C"]

            self.getCellBarcode = self._getCellBarcodeString
            self.getBarcodes = self._getBarcodesString

    def getReadCounts(self):
        return self.read_counts

    def __call__(self, read1, read2=None):

        self.read_counts['Input Reads'] += 1

        umi_values = self.getBarcodes(read1, read2)

        if umi_values is None:
            return None
        else:
            cell, umi, new_seq, new_seq2 = umi_values

        #if self.filter_sgrna:
        cell, gene = self.filterCellBarcode(cell)

        if cell is None:
                return None

        self.read_counts['Reads output'] += 1

        # if UMI could be on either read, use umi_values to identify
        # which read(s) it was on
        # Otherwise, use input from user to identiy which reads need updating

        new_identifier = addBarcodesToIdentifier(read1, umi, cell)
        read1.identifier = new_identifier
        if self.pattern:  # seq and quals need to be updated
            read1.seq = new_seq

        if read2:
            new_identifier2 = addBarcodesToIdentifier(read2, umi, cell)
            read2.identifier = new_identifier2
            if self.pattern2:   # seq and quals need to be updated
                read2.seq = new_seq2

        if read2 is None:
            return gene, read1
        else:
            return gene, read1, read2

class UMIClusterer:
    '''A functor that clusters a dictionary of UMIs and their counts.
    The primary return value is either a list of representative UMIs
    or a list of lists where each inner list represents the contents of
    one cluster.

    Optionally:

      - identify the parent UMIs and return:
         - selected reads
         - umis
         - counts

    The initiation of the functor defines the methods:

      ** get_adj_list ** - returns the edges connecting the UMIs

      ** get_connected_components ** - returns clusters of connected components
                                       using the edges in the adjacency list

      ** get_groups ** - returns the groups of umis,
                         with the parent umi at position 0

    Note: The get_adj_list and connected_components methods are not required by
    all custering methods. Where there are not required, the methods return
    None or the input parameters.

    '''

    def _get_adj_list_directional(self, umis, counts, threshold):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1.

        counts is a direction with keys are umis and values are counts of that umi
        '''
        adj_list = {umi: [] for umi in umis}
        if len(umis) > 25:
            umi_length = len(umis[0])
            substr_idx = build_substr_idx(umis, umi_length, threshold)
            iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
        else:
            iter_umi_pairs = itertools.combinations(umis, 2)
        for umi1, umi2 in iter_umi_pairs:
            if edit_distance(umi1, umi2) <= threshold:
                if counts[umi1] >= (counts[umi2] * 2) - 1:
                    adj_list[umi1].append(umi2)
                if counts[umi2] >= (counts[umi1] * 2) - 1:
                    adj_list[umi2].append(umi1)

        return adj_list

    def _get_adj_list_null(self, umis, counts, threshold):
        ''' for methods which don't use a adjacency dictionary'''
        return None

    def _get_connected_components_adjacency(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        # TS: TO DO: Work out why recursive function doesn't lead to same
        # final output. Then uncomment below

        # if len(graph) < 10000:
        #    self.search = breadth_first_search_recursive
        # else:
        #    self.search = breadth_first_search

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                # component = self.search(node, graph)
                component = breadth_first_search(node, graph)
                found.update(component)
                components.append(component)
        return components

    def _get_connected_components_null(self, umis, adj_list, counts):
        ''' for methods which don't use a adjacency dictionary'''
        return umis

    # "group" methods #

    def _group_directional(self, clusters, adj_list, counts):
        ''' return groups for directional method'''

        observed = set()
        groups = []
        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))
                observed.update(cluster)
            else:
                cluster = sorted(cluster, key=lambda x: counts[x],
                                 reverse=True)
                # need to remove any node which has already been observed
                temp_cluster = []
                for node in cluster:
                    if node not in observed:
                        temp_cluster.append(node)
                        observed.add(node)
                groups.append(temp_cluster)

        return groups

    def _group_cluster(self, clusters, adj_list, counts):
        ''' return groups for cluster or directional methods'''

        groups = []
        for cluster in clusters:
            groups.append(sorted(cluster, key=lambda x: counts[x],
                                 reverse=True))

        return groups

    def __init__(self, cluster_method="directional"):
        ''' select the required class methods for the cluster_method'''

        self.max_umis_per_position = 0
        self.total_umis_per_position = 0
        self.positions = 0

        if cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_directional

    def __call__(self, umis, threshold):
        '''umis is a dictionary that maps UMIs to their counts'''

        counts = umis
        umis = list(umis.keys())

        self.positions += 1

        number_of_umis = len(umis)

        self.total_umis_per_position += number_of_umis

        if number_of_umis > self.max_umis_per_position:
            self.max_umis_per_position = number_of_umis

        len_umis = [len(x) for x in umis]

        assert max(len_umis) == min(len_umis), (
                "not all umis are the same length(!):  %d - %d" % (
            min(len_umis), max(len_umis)))

        adj_list = self.get_adj_list(umis, counts, threshold)
        clusters = self.get_connected_components(umis, adj_list, counts)
        final_umis = [list(x) for x in
                      self.get_groups(clusters, adj_list, counts)]

        return final_umis

def fastqIterate(infile):
    '''iterate over contents of fastq file.'''

    def convert2string(b):
        if type(b) == str:
            return b
        else:
            return b.decode("utf-8")

    while 1:
        line1 = convert2string(infile.readline())
        if not line1:
            break
        if not line1.startswith('@'):
            error("parsing error: expected '@' in line %s" % line1)
        line2 = convert2string(infile.readline())
        line3 = convert2string(infile.readline())
        if not line3.startswith('+'):
            error("parsing error: expected '+' in line %s" % line3)
        line4 = convert2string(infile.readline())
        # incomplete entry
        if not line4:
            error("incomplete entry for %s" % line1)

        yield Record(line1[1:-1], line2[:-1], line4[:-1])

def joinedFastqIterate(fastq_iterator1, fastq_iterator2, strict=True):
    '''This will return an iterator that returns tuples of fastq records.
    At each step it will confirm that the first field of the read name
    (before the first whitespace character) is identical between the
    two reads. The response if it is not depends on the value of
    :param:`strict`. If strict is true an error is returned. If strict
    is `False` the second file is advanced until a read that matches
    is found.
    This allows for protocols where read one contains cell barcodes, and these
    reads have been filtered and corrected before processing without regard
    to read2
    '''

    for read1 in fastq_iterator1:
        read2 = next(fastq_iterator2)
        pair_id = read1.identifier.split()[0]
        if not strict:
            while read2.identifier.split()[0] != pair_id:
                read2 = next(fastq_iterator2)
        if not read2.identifier.split()[0] == pair_id:
            raise ValueError("\nRead pairs do not match\n%s != %s" %
                             (pair_id, read2.identifier.split()[0]))
        yield (read1, read2)

def get_cell_umi_read_string(read_id, sep="_"):
    ''' extract the umi and cell barcode from the read id (input as a
    string) using the specified separator '''

    try:
        return (read_id.split(sep)[-1].encode('utf-8'),
                read_id.split(sep)[-2].encode('utf-8'))
    except IndexError:
        raise ValueError(
            "Could not extract UMI or CB from the read ID, please"
            "check UMI and CB are encoded in the read name:"
            "%s" % read_id)

def get_gene_count_tab(df_records, bc_getter=None):
    gene = None
    counts = collections.Counter()
    df_records.sort_values(axis=0, by="gene", inplace=True)

    for row in df_records.itertuples(index=False):
        read_id, assigned_gene = list(row)

        if assigned_gene != gene:
            if gene:
                yield gene, counts

            gene = assigned_gene
            counts = collections.defaultdict(collections.Counter)

        cell, umi = bc_getter(read_id)
        counts[cell][umi] += 1

    yield gene, counts

def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched

def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    args.sgrna = argv[0]
    args.barcode = argv[1]
    args.whitelist = argv[2]
    args.output = argv[3]

    print("=====================")
    print(args)
    print("======================")

    #parser.add_argument("-s", "--sgrna", dest="sgrna", help="Fastq file contains sgRNAs")
    #parser.add_argument("-b", "--barcode", dest="barcode", help="Fastq file contains barcode")
    #parser.add_argument("-w", "--whitelist", dest="whitelist", help="Valid sgRNA file in CSV format: id,sgRNA,gene")
    #parser.add_argument("-t", "--threshold", dest="threshold", help="Edit-distance for barcode collapsing. Default 1.", type=int, default=1)
    #parser.add_argument("-o", "--output", dest="output", help="Output file name", default="count_umi_output.txt")

    if len(argv) == 1:
        parser.print_help()
        sys.exit()



    fout = args.output
    #if os.path.isfile(fout):
    #    overwrite = input('\nFile already exists. Overwrite? Y/N ')
    #    if overwrite.lower() == 'n':
    #        exit(0)

    read_out = openFile(fout, "w")

    if not args.whitelist:
        error("must provide a whitelist (--whitelist)")

    read1s = fastqIterate(openFile(args.sgrna,'r'))
    read2s = fastqIterate(openFile(args.barcode,'r'))

    ReadExtractor = ExtractFilterAndUpdate(method="string",
                                           pattern="CCCCCCCCCCCCCCCCCCCC",
                                           pattern2="NNNNNNNNNNNNNNNNNNNN",
                                           extract_cell=True,
                                           filter_sgrna=True
                                           )

    cell_whitelist, gene_whitelist = getUserDefinedBarcodes(args.whitelist)
    ReadExtractor.cell_whitelist = cell_whitelist
    ReadExtractor.gene_whitelist = gene_whitelist

    records = []

    for read1, read2 in joinedFastqIterate(read1s, read2s):
        reads = ReadExtractor(read1, read2) ## type(reads): tuple
        if not reads:
            continue
        else:
            gene, new_read1, new_read2 = reads

        records.append((new_read1.identifier.split(" ")[0], gene))

    records = pd.DataFrame(records, columns =['read_id','gene'])

    # start counting
    nInput, nOutput = 0, 0
    bc_getter = partial(get_cell_umi_read_string, sep="_")
    processor = UMIClusterer("directional")
    for gene, counts in get_gene_count_tab(
            records,
            bc_getter=bc_getter):

        for cell in counts.keys():

            umis = counts[cell].keys()

            nInput += sum(counts[cell].values())

            # group the umis
            groups = processor(counts[cell], threshold=args.threshold)
            gene_count = len(groups)
            for group in groups:
                group.sort()
                umi_ids = ",".join(u.decode("utf-8") for u in group)
                umi_seq= str(group[0].decode("utf-8"))
                umi_counts = 0
                for umi in group:
                    umi_counts += counts[cell].get(umi)

                read_out.write("%s\t%s\t%i\n" %(cell.decode("utf-8")+"_"+umi_seq, gene+"_"+cell.decode("utf-8"), umi_counts))

            nOutput += gene_count

    logging.info("Number of reads counted: %i" % nOutput)


if __name__ == "__main__":
    argv = [snakemake.input[0], snakemake.input[1], snakemake.input[2], snakemake.output[0]]
    sys.exit(main(argv))
