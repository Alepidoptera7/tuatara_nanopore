import subprocess
from Bio.Align.Applications import ClustalOmegaCommandline

class read_to_gene:
    """
    This class develops genomic categorization.
    """

    def __init__(self):
        self.gene_file_list = ["LadyAliceND4.fasta", "LadyAliceND4.fasta", "LadyAliceND5.fasta", "LadyAliceCO1.fasta", "LadyAliceCO2.fasta",
                               "StephensIslandCO1.fasta", "StephensIslandCO2.fasta", "StephensIslandND4.fasta", "StephensIslandND5.fasta"]
        self.file_dict = {i: [] for i in self.gene_file_list}
        self.percent_id_dict = {}

    def nanopore_parser(self):
        """This def is designed to open a given file and distribute genomic sequences to the k-means algorithm.

        Input:read file
        Output: read file data
        """

        nanopore_string = ""
        seq_num = 0
        with open("SLZ14846_nanopore_1k-18k.fasta") as in_file:
            for i in range(2000):
                line = next(in_file)
                nanopore_string += line

                if line == "\n":
                    split_nanopore_string = nanopore_string.split("@v3.5.2")
                    header, seq = split_nanopore_string[0], split_nanopore_string[1]
                    nanopore_string = ""

                    if header and seq != "\n" or "":
                        self.alignment_file_creator(header, seq)
                        seq_num += 1
                        print("seq num: ", seq_num)

            in_file.close()


    def alignment_file_creator(self, header, seq):
        """"""

        for ref_file in self.gene_file_list:
            with open(ref_file) as ref:
                ref_seq = ""
                ref_header = ""
                file_to_aln = open("clustal_infile.fa", "w")

                for line in ref:
                    if ">" in line:
                        ref_header = line.strip("\n")
                    else:
                        ref_seq += line

                if len(ref_seq) > len(seq):
                    file_to_aln.write(header)
                    file_to_aln.write("\n")
                    file_to_aln.write("$" + seq + "$")
                    file_to_aln.write("\n")
                    file_to_aln.write(ref_header + " " + ref_file)
                    file_to_aln.write("\n")
                    file_to_aln.write(ref_seq)

                if len(ref_seq) < len(seq):
                    file_to_aln.write(ref_header + " " + ref_file)
                    file_to_aln.write("\n")
                    file_to_aln.write("$" + ref_seq + "$")
                    file_to_aln.write("\n")
                    file_to_aln.write(header)
                    file_to_aln.write("\n")
                    file_to_aln.write(seq)

            cline = self.biopython_clustalw("clustal_infile.fa")
            self.sub_process(cline)
            print(header)
            print(ref_header + " " + ref_file)
            percent_id = self.percent_id_calculator()

            if percent_id > 90:
                self.file_dict[ref_file].append((header, percent_id))

            ref.close()

    def biopython_clustalw(self, infile):
        """The purpose of this def is to develop a command to call clustal command line tool."""

        clustalOmega_exe = r"C:/Users/Quin The Conquoror!/Desktop/clustal-omega-1.2.2-win64/clustalo"
        cline_outfile = "cline_outfile.fa"
        cline = ClustalOmegaCommandline(clustalOmega_exe, infile=infile, outfile=cline_outfile, outfmt="fasta", verbose=True, force=True)
        cline_str = str(cline)

        return cline_str

    def sub_process(self, cline):
        subprocess.run(cline)

    def percent_id_calculator(self):
        """The purpose of this def is to caculate percent identity between
         a given sequence and the reference genome from previously generated alignments.

         Input: previously generated .fasta alignment file
         Output: percent identy per alignment
         """

        alignment_string = ""
        alignment_string_list = []
        gap_count = 0

        with open("cline_outfile.fa") as aln:
            for line in aln:
                alignment_string_list.append(line)
                #print(line)

        for line in alignment_string_list[1:]:
            alignment_string += line

        start_char_index = alignment_string.index("$")
        end_char_index = alignment_string.rindex("$")

        for base in alignment_string[start_char_index:end_char_index]:
            if base == "-":
                gap_count += 1

        aln_len = end_char_index - start_char_index
        match_count = aln_len - gap_count
        percent_identity = ((match_count * 100) / aln_len)

        if percent_identity > 98:
            print(alignment_string_list[0])
            print(alignment_string[start_char_index:end_char_index])
            print("start_char_index", start_char_index)
            print("end char index", end_char_index)
            #print("gap count", gap_count)
            #print("alignment length", aln_len)
            #print("match count", match_count)
            print("Percent Identity: ", ((match_count * 100)/aln_len))

        return percent_identity

    def output_file_creator(self):
        """This def is designed to take the high percent identiy reads for each gene alignment and write them into separate files for review."""
        for key in self.file_dict.keys():
            print(key, self.file_dict[key])
            new_file = open(key + "output.fa", "w")
            self.file_dict[key].sort(key=lambda a: a[1])
            for result in self.file_dict[key]:
                    new_file.write(str(result))
                    new_file.write("\n")

    def driver(self):
        """Driver calls functions."""
        self.nanopore_parser()
        self.output_file_creator()


def main():
    class_access = read_to_gene()
    class_access.driver()

if __name__ == '__main__':
    main()