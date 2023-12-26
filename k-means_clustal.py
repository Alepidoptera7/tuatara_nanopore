import random
import math
import matplotlib.pyplot as plt
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import subprocess
from Bio.Align.Applications import ClustalOmegaCommandline

class k_means:
    """
    This class develops genomic categorization.
    """

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        self.read_dict = {}
        self.k = 10
        self.kmers = []

        self.multi_read_cluster_holder_dict = {}
        self.centroids_dict = {}
        self.multi_read_centroid_change_list = [i for i in range(self.k)]
        self.seq_dict = {}
        self.sequence_assembly_dict = {}

        self.x_coords_dict = {}
        self.y_coords_dict = {}

        self.random_string = ''
        self.nanopore_read_dict = {}
        self.percent_id_list = []
        self.percent_id_dict = {}

    def nanopore_parser(self):
        """This def is designed to open a given file and distribute genomic sequences to the k-means algorithm.

        Input:read file
        Output: read file data
        """
        read_count = 0
        nanopore_string = ""

        with open("SLZ14846-nanopore-test.fastq") as in_file:
            for line in in_file:

                if "@" and "read" in line:
                    nanopore_string += line

                if "@" not in line:
                    nanopore_string += line

                if line == "+\n":
                    split_nanopore_string = nanopore_string.split("\n")
                    header, seq = split_nanopore_string[0], split_nanopore_string[1]
                    read_count += 1
                    print(read_count, header)
                    self.pr_matrix(header, seq)
                    nanopore_string = ""

    def pr_matrix(self, head, seq):
        """

        Develops a dictionary in which each kmer key is associated with a list of unique (x,y) = (position, frequency) coordinate tuples.

        These values are attached to the given kmer and used at y-coordinates for Euclidian distance calculations.

        1) The sequence is parsed into 4-mers.
        2) The 4-mers are used for keys in the dictionary.
        3) The number of times each 4-mer occurs in the sequence is counted and divided by the number of possible 4mers, 256.


        Input: FASTA sequence
        Output: Dictionary with frequency Pr values attached to kmers.
        """

        self.kmers = [seq[i:i + 4] for i in range(0, len(seq), 4)]
        self.seq_dict[head] = seq
        self.sequence_assembly_dict[head] = {}

        kmer_xy_coordinate_list = []
        for i in range(0, len(self.kmers)):
            if len(self.kmers[i]) == 4:
                kmer_xy_coordinate_list.append((i*4, self.kmers.count(self.kmers[i]) / 256))

        self.read_dict[head] = kmer_xy_coordinate_list

    def random_centroids(self):
        """Randomly selects k number of initial centroids by position in FASTA sequence.

        The initial centroids are attached to the read in the overarching read_dict.

        These centroids will be updated and used for cluster identification.

        Input: none
        Output: dictionary random centroids

        """

        for read in self.read_dict.keys():
            random_centroids = []

            for i in range(0, self.k):
                random_centroids.append(self.read_dict[read][random.randint(0, self.k)])

            self.centroids_dict[read] = random_centroids
            self.multi_read_cluster_holder_dict[read] = {i: [] for i in range(0, self.k)}

    def distance_calculation(self):
        """Euclidian distance is calculated between each centroid and all non-self points.

        Input: a list of points and centroids.
        Output: a list of distances between each centroid and all other points.
        """

        for read in self.read_dict.keys():
            temp_placement_dict = {i: [] for i in range(0, self.k)}

            for point in self.read_dict[read]:
                distance_list = []
                x2, y2 = point[0], point[1]

                for centroid in self.centroids_dict[read]:
                    x1, y1 = centroid[0], centroid[1]

                    #distance calculation between centroid and point
                    distance = math.sqrt(abs(x2 - x1) + abs(y2 - y1))
                    distance_list.append(distance)

                temp_placement_dict[distance_list.index(min(distance_list))].append(point)

            for i in range(self.k):
                self.multi_read_cluster_holder_dict[read][i] = temp_placement_dict[i]

    def centroid_recalculation(self):
        """Recalculates centroids by cluster averaging.

        Input: Current cluster information
        Output: List of adjusted centroids as component to self.read_dict[read][1]
        """

        for read in self.multi_read_cluster_holder_dict.keys():
            temp_centroid_list = []
            centroid_change_list = []
            for i in range(self.k):
                if len(self.multi_read_cluster_holder_dict[read][i]) > 0:
                    new_centroid_x = sum([tup[0] for tup in self.multi_read_cluster_holder_dict[read][i]]) / len(self.multi_read_cluster_holder_dict[read][i])
                    new_centroid_y = sum([tup[1] for tup in self.multi_read_cluster_holder_dict[read][i]]) / len(self.multi_read_cluster_holder_dict[read][i])

                    temp_centroid_list.append((new_centroid_x, new_centroid_y))
                    centroid_change_list.append(abs(new_centroid_x - self.centroids_dict[read][i][0]) + abs(new_centroid_y - self.centroids_dict[read][i][1]))

            self.multi_read_centroid_change_list = centroid_change_list
            self.centroids_dict[read] = temp_centroid_list

    def sequence_assembly(self):
        """Assembling the sequences represented by the predicted clusters.

        Input: Finished Clusters
        Output: Separated genomic sequences
        """

        for header in self.multi_read_cluster_holder_dict.keys():
            print(header)

            self.x_coords_dict[header] = {i: [] for i in range(0, self.k)}
            self.y_coords_dict[header] = {i: [] for i in range(0, self.k)}

            for i in self.multi_read_cluster_holder_dict[header].keys():
                cluster_sequence = ""
                self.sequence_assembly_dict[header][i] = ""

                x_coords_list = []
                y_coords_list = []

                for point in self.multi_read_cluster_holder_dict[header][i]:
                    cluster_sequence += self.seq_dict[header][point[0]: point[0] + 4]

                    x_coords_list.append(point[0])
                    y_coords_list.append(point[1])

                self.sequence_assembly_dict[header][i] = cluster_sequence

                self.x_coords_dict[header][i] = x_coords_list
                self.y_coords_dict[header][i] = y_coords_list

                if cluster_sequence != "":

                    print("cluster: ", i, " cluster size: ", "(", len(self.multi_read_cluster_holder_dict[header][i]), "/", len(self.read_dict[header]), ")",
                          " portion of read: ",
                         len(self.multi_read_cluster_holder_dict[header][i]) / len(self.read_dict[header]))

                    print(self.sequence_assembly_dict[header][i])

                    #self.blast(cluster_sequence)
                    self.seq_to_fasta(header, cluster_sequence)

                    print("_____________________________________________________________________________\n")
                    print("_____________________________________________________________________________\n")


    def seq_to_fasta(self, header, seq):
        """The purpose of this def is to create a fasta file from the provided sequence from the original file.
        The same file name will be used so that a multitude of files will not be maintained.
        """

        infile = open("current_aln.fa", "w")
        infile.write(header)
        infile.write(seq)
        self.biopython_clustalw(infile)

    def biopython_clustalw(self, infile):
        """The purpose of this def is to develop a command to call clustal command line tool."""

        clustalOmega_exe = r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        cline_outfile = "current_aln.fa"
        print(cline_outfile)

        cline = ClustalOmegaCommandline(clustalOmega_exe, infile=infile, verbose=True, outfile=cline_outfile, outfmt="fasta")
        cline_str = str(cline)
        subprocess.run(cline_str)

        self.percentid_calculator(cline_outfile)

    def percentid_calculator(self, cline_outfile):
        """The purpose of this def is to caculate percent identity between
         a given sequence and the reference genome from previously generated alignments.

         Input: previously generated .fasta alignment file
         Output: percent identy per alignment
         """

        base_list = ["N", "A", "G", "C", "T", "n"]
        alignment_string = ""
        base_index_list = []
        alignment_string_list = []
        gap_count = 0

        with open(cline_outfile) as aln:
            for line in aln:
                alignment_string_list.append(line)

        for line in alignment_string_list[1:]:
            alignment_string += line

        for base in base_list:
            base_index_list.append(alignment_string.index(base))

        base_index_list.sort()

        for base in alignment_string[base_index_list[0]:base_index_list[-1]]:
            if base == "-":
                gap_count += 1

        aln_len = base_index_list[-1] - base_index_list[0]
        match_count = aln_len - gap_count

        self.percent_id_dict[cline_outfile] = ((match_count * 100) / aln_len)

        percent_identity = ((match_count * 100) / aln_len)

        if percent_identity >= 0.95:
            self.percent_id_list.append((cline_outfile, percent_identity))

        # print(aln_outfile)
        # print("base index list", base_index_list)
        # print("gap count", gap_count)
        # print("alignment length", aln_len)
        # print("match count", match_count)
        # print("Percent Identity: ", ((match_count * 100)/aln_len))


    def blast(self, cluster_sequence):
        """Develops a BLAST request for each clustered sequence via BioPython.

        Design taken from Biopython Manual.
        """
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence=cluster_sequence, perc_ident=85)

        blast_record = NCBIXML.read(result_handle)

        count = 0

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                count += 1
                print(f"sequence: {alignment.title}, length: {alignment.length}nt, e-value: {hsp.expect}")
                print('Percent Identity: %', round(hsp.identities / hsp.align_length * 100, 2))

        print(f"There are {count} similar sequences in Blast output\n")

    def driver(self):
        """Drivers iterative processes"""

        self.nanopore_parser()
        self.random_centroids()

        self.distance_calculation()

        iteration = 0
        while sum(self.multi_read_centroid_change_list) > 0:
            self.centroid_recalculation()
            self.distance_calculation()

            iteration += 1

        print("total iterations of k-means clustering: ", iteration)

        self.sequence_assembly()


def main():
    class_access = k_means()
    class_access.driver()

if __name__ == '__main__':
    main()