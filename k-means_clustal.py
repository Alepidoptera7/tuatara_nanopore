import random
import math
import numpy
import subprocess
from Bio.Align.Applications import ClustalOmegaCommandline

class k_means:
    """
    This class develops genomic categorization.
    """

    def __init__(self):
        '''contructor: saves attribute fname '''
        self.read_dict = {}
        self.k = 1
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

        Develops a dictionary in which each kmer key is associated with a list of unique (x,y) = (position, frequency) coordinate tuples for 3 reading frames.

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
                x_coords_list, y_coords_list = [], []

                for point in self.multi_read_cluster_holder_dict[header][i]:
                    cluster_sequence += self.seq_dict[header][point[0]: point[0] + 4]
                    x_coords_list.append(point[0])
                    y_coords_list.append(point[1])

                self.x_coords_dict[header][i], self.y_coords_dict[header][i]= x_coords_list, y_coords_list

                if cluster_sequence != "":
                    print("cluster: ", i, " cluster size: ", "(", len(cluster_sequence), "/", len(self.seq_dict[header]), ")",
                          " portion of read: ", len(self.multi_read_cluster_holder_dict[header][i]) / len(self.read_dict[header]))

                    print(cluster_sequence)

                    cluster_sequence = cluster_sequence + "\n"
                    self.seq_to_fasta(header, cluster_sequence)

                    print("_____________________________________________________________________________\n")

    def seq_to_fasta(self, header, seq):
        """The purpose of this def is to create a fasta file from the provided sequence from the original file.
        The same file name will be used so that a multitude of files will not be maintained.
        """
        print("header, seq:")
        print(header, seq)
        new_file = open("infile.fa", "w")
        new_file.write(">")
        new_file.write(header)
        new_file.write("\n")
        new_file.write(seq)
        new_file.write("\n")

        with open("MN864230.1.fa") as ref:
            for line in ref:
                new_file.write(line)

        new_file.close()

        cline = self.biopython_clustalw("infile.fa")
        self.sub_process(cline)
        self.percentid_calculator(header, "cline_aln", len(seq))

    def biopython_clustalw(self, infile):
        """The purpose of this def is to develop a command to call clustal command line tool."""

        clustalOmega_exe = r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        cline_outfile = "cline_aln"
        cline = ClustalOmegaCommandline(clustalOmega_exe, infile=infile, outfile=cline_outfile, outfmt="fasta", verbose=True, force=True )
        cline_str = str(cline)

        return cline_str

    def sub_process(self, cline):
        subprocess.run(cline)

    def percentid_calculator(self, header, cline_outfile, seq_len):
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

        #print("base index list", base_index_list)
        #print("gap count", gap_count)
        #print("alignment length", aln_len)
        #print("match count", match_count)
        print("Percent Identity: ", ((match_count * 100)/aln_len))

        self.percent_id_list.append((header, percent_identity, seq_len))

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

        self.percent_id_list.sort(key=lambda a: a[1])

        for i in self.percent_id_list:
            print(i)


def main():
    class_access = k_means()
    class_access.driver()

if __name__ == '__main__':
    main()