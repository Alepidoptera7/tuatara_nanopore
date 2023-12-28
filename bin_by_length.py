class bin_by_len:
    """
    This class develops genomic categorization.
    """

    def __init__(self):
        self.read_bin_dict = {i:[] for i in range(0, 200000, 1000)}

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
                    #header_split = header.split("read=")
                    #print(header_split)
                    read_count += 1
                    self.bin_placement(header, seq)
                    nanopore_string = ""

    def bin_placement(self, header, seq):
        """This def places the reads in the bin corresponding to their length. """

        seq_len = len(seq)

        for i in self.read_bin_dict.keys():
            if i - 1000 < seq_len < i:
                #print(i, header, seq_len)
                self.read_bin_dict[i].append((header, seq_len))

    def print_out(self):
        """Prints the bin dict."""

        for i in self.read_bin_dict.keys():
            print(len(self.read_bin_dict[i]), self.read_bin_dict[i])

    def driver(self):
        """Drivers iterative processes"""

        self.nanopore_parser()
        self.print_out()

def main():
    class_access = bin_by_len()
    class_access.driver()

if __name__ == '__main__':
    main()