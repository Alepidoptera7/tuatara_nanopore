import subprocess
from Bio.Align.Applications import ClustalOmegaCommandline

class nanopore_clustal:

    def __init__(self):
        self.nanopore_read_dict = {}
        self.aln_outfile_list = []
        self.percent_id_list = []
        self.percent_id_dict = {}

    def nanopore_parser(self):
        """This def is designed to open a given file and distribute genomic sequences to the k-means algorithm.

        Input:read file
        Output: read file data
        """
        read_count = 0
        nanopore_string = ""
        with open("SLZ14846-nanopore.fastq") as in_file:
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

                    if header not in self.nanopore_read_dict.keys():
                        self.nanopore_read_dict[header] = seq
                    else:
                        self.nanopore_read_dict[header] += seq

                    nanopore_string = ""

    def seq_to_fasta(self):
        """The purpose of this def is to create a fasta file from the provided sequence from the original file.
        The same file name will be used so that a multitude of files will not be maintained.
        """

        for header in self.nanopore_read_dict.keys():
            infile = open("current_aln.fa", "w")
            infile.write(header)
            infile.write(self.nanopore_read_dict[header])
            self.biopython_clustalw(infile)

    def biopython_clustalw(self, infile):
        """The purpose of this def is to develop a command to call clustal command line tool."""

        clustalOmega_exe = r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        cline_outfile = "current_aln"
        print(cline_outfile)

        cline = ClustalOmegaCommandline(clustalOmega_exe, infile=infile, verbose=True, outfile=cline_outfile, outfmt="fasta", percentid=True)

        return str(cline)

    def sub_process(self, cline):
        subprocess.run(cline)

    def percentid_calculator(self):
        """The purpose of this def is to caculate percent identity between
         a given sequence and the reference genome from previously generated alignments.

         Input: previously generated .fasta alignment file
         Output: percent identy per alignment
         """

        base_list = ["N", "A", "G", "C", "T", "n"]

        for aln_outfile in self.aln_outfile_list:
            alignment_string = ""
            base_index_list = []
            alignment_string_list = []
            gap_count = 0
            with open(aln_outfile) as aln:
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

            self.percent_id_dict[aln_outfile] = ((match_count * 100) / aln_len)

            percent_identity = ((match_count * 100) / aln_len)

            if percent_identity >= 0.95:

                self.percent_id_list.append((aln_outfile, percent_identity))

            # print(aln_outfile)
            # print("base index list", base_index_list)
            # print("gap count", gap_count)
            # print("alignment length", aln_len)
            # print("match count", match_count)
            # print("Percent Identity: ", ((match_count * 100)/aln_len))

    def driver(self):
        self.nanopore_parser()
        self.percentid_calculator()

        self.percent_id_list.sort(key=lambda a: a[1])

        for i in self.percent_id_list:
            print(i)


def main():
    class_access = nanopore_clustal()
    class_access.driver()

if __name__ == "__main__":
    main()