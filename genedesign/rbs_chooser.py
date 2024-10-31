from typing import Set
from genedesign.models.rbs_option import RBSOption
import pandas as pd # type: ignore

from genedesign.seq_utils import Translate, hairpin_counter
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance 

pruned_data = pd.read_csv('pruned_data.csv')
genes_info = pd.read_csv('genes_info.csv')


class RBSChooser:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """


    # Define class variables
    rbs_options: Set[RBSOption]
    translator: Translate

    def initiate(self) -> None:
        """
        Initialization method for RBSChooser.
        """
        # Instantiate the variables
        self.rbs_options = []
        self.translator = Translate()
        self.translator.initiate()

        # Create RBSOption instances
        for locus_tag, abundance in pruned_data:
            if locus_tag in genes_info:
                gene_info = genes_info[locus_tag]
                utr = gene_info['UTR']
                cds = gene_info['CDS']
                gene_name = gene_info['gene']
                first_six_aas = self.translator.run(cds[:18])

                rbs_option = RBSOption(utr=utr, cds=cds, gene_name=gene_name, first_six_aas=first_six_aas)
                self.rbs_options.append(rbs_option)

    def score_rbs_option(self, rbs_option: RBSOption, input_cds: str) -> float:
        """
        Score an RBSOption based on its compatibility with the input CDS.

        Parameters:
            rbs_option (RBSOption): The RBSOption to score.
            input_cds (str): The input coding sequence.

        Returns:
            float: A score representing the compatibility of the RBS option with the input CDS (lower is better).
        """
        # Calculate edit distance
        input_first_six_aas = self.translator.run(input_cds[:18])
        edit_distance = calculate_edit_distance(input_first_six_aas, rbs_option.first_six_aas)


        # Calculate hairpin count
        combined_sequence = rbs_option.utr + input_cds[:50]  # Consider first 50 nucleotides of input CDS
        hairpin_count = hairpin_counter(combined_sequence)

        # Calculate score (lower is better)
        score = edit_distance + hairpin_count * 2 # Weigh hairpin count more heavily

        return score


    def filter_ignored_options(self, ignores: Set[RBSOption]) -> Set[RBSOption]:
        """
        Filters out RBSOptions that are in the ignores set.

        Parameters:
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore.

        Returns:
        - List[RBSOption]: A list of RBSOption instances not in the ignores set.
        """

        return [option for option in self.rbs_options if option not in ignores]


    def run(self, cds: str, ignores: Set[RBSOption]) -> RBSOption:
        """
        Executes the RBS selection process for the given CDS.

        Parameters:
        - cds (str): The coding sequence to pair with an RBS.
        - ignores (Set[RBSOption]): A set of RBSOption instances to ignore during selection.

        Returns:
        - RBSOption: The selected RBSOption that best pairs with the given CDS.
        """

        filtered_options = self.filter_ignored_options(ignores)

        if not filtered_options:
            return None

        # Score each RBS option
        scored_options = [(option, self.score_rbs_option(option, cds)) for option in filtered_options]

        # Sort by score (lower is better)
        scored_options.sort(key=lambda x: x[1])

        # Return the best-scoring option
        return scored_options[0][0]

    
if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT..."

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)



'''

class RBSChooser:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """

    def __init__(self):
        self.rbsOptions = []

    def initiate(self) -> None:
        """
        Populates the RBS options list with predefined options.
        """
        # Add RBS options based on their sequences and properties
        opt1 = RBSOption(
            utr="aaagaggagaaatactag",
            cds="atggcttcctccgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaa",
            gene_name="BBa_b0034",
            first_six_aas="MASSED"
        )
        opt2 = RBSOption(
            utr="tcacacaggaaagtactag",
            cds="atgactcaacgtatcgcatatgtaactggtggtatgggtggtatcggtactgcaatttgccagcgtctggcgaaagacggtttccgtgttgttgcgggctgcggtccgaactccccgcgtcgtgaaaagtggctggaacaacagaaagccctgggcttcgacttcattgcctccgagggtaatgtagctgactgggattccaccaagactgccttcgataaagttaaatctgaagtgggcgaagtagatgtactgatcaacaacgccggtattactcgtgatgtcgtattccgcaaaatgacccgtgcagactgggatgcagttatcgacaccaacctgacgtctctgttcaacgttaccaaacaggttattgatggtatggctgaccgtggctggggccgcatcgtgaacatctctagcgttaacggccaaaaaggccaatttggtcagacgaattacagcacggctaaagcaggcctgcacggtttcaccatggcactggcgcaggaagtggcgaccaaaggtgttaccgttaataccgtttctccaggttacatcgccaccgatatggttaaggctatccgccaagatgttctggacaagatcgtggctaccattccggttaaacgcctgggcctgccggaagaaattgcgtccatctgtgcgtggctgagctccgaagagtctggtttttccaccggtgcggatttctctctgaacggtggtctgcacatgggttga",
            gene_name="BBa_b0032",
            first_six_aas="MTQRIA"
        )
        opt3 = RBSOption(
            utr="CCATACCCGTTTTTTTGGGCTAACAGGAGGAATTAAcc",
            cds="atgGacacAattaacatcgctaagaacgacttctctgacatcgaactggctgctatcccgttcaacactctggctgaccattacggtgagcgtttagctcgcgaacagttggcccttgagcatgagtcttacgagatgggtgaagcacgcttccgcaagatgtttgagcgtcaacttaaagctggtgaggttgcggataacgctgccgccaagcctctcatcactaccctactccctaagatgattgcacgcatcaacgactggtttgaggaagtgaaagctaagcgcggcaagcgcccgacagccttccagttcctgcaagaaatcaagccggaagccgtagcgtacatcaccattaagaccactctggcttgcctaaccagtgctgacaatacaaccgttcaggctgtagcaagcgcaatcggtcgggccattgaggacgaggctcgcttcggtcgtatccgtgaccttgaagctaagcacttcaagaaaaacgttgaggaacaactcaacaagcgcgtagggcacgtctacaagaaagcatttatgcaagttgtcgaggctgacatgctctctaagggtctactcggtggcgaggcgtggtcttcgtggcataaggaagactctattcatgtaggagtacgctgcatcgagatgctcattgagtcaaccggaatggttagcttacaccgccaaaatgctggcgtagtaggtcaagactctgagactatcgaactcgcacctgaatacgctgaggctatcgcaacccgtgcaggtgcgctggctggcatctctccgatgttccaaccttgcgtagttcctcctaagccgtggactggcattactggtggtggctattgggctaacggtcgtcgtcctctggcgctggtgcgtactcacagtaagaaagcactgatgcgctacgaagacgtttacatgcctgaggtgtacaaagcgattaacattgcgcaaaacaccgcatggaaaatcaacaagaaagtcctagcggtcgccaacgtaatcaccaagtggaagcattgtccggtcgaggacatccctgcgattgagcgtgaagaactcccgatgaaaccggaagacatcgacatgaatcctgaggctctcaccgcgtggaaacgtgctgccgctgctgtgtaccgcaaggacaaggctcgcaagtctcgccgtatcagccttgagttcatgcttgagcaagccaataagtttgctaaccataaggccatctggttcccttacaacatggactggcgcggtcgtgtttacgctgtgtcaatgttcaacccgcaaggtaacgatatgaccaaaggactgcttacgctggcgaaaggtaaaccaatcggtaaggaaggttactactggctgaaaatccacggtgcaaactgtgcgggtgtcgataaggttccgttccctgagcgcatcaagttcattgaggaaaaccacgagaacatcatggcttgcgctaagtctccactggagaacacttggtgggctgagcaagattctccgttctgcttccttgcgttctgctttgagtacgctggggtacagcaccacggcctgagctataactgctcccttccgctggcgtttgacgggtcttgctctggcatccagcacttctccgcgatgctccgagatgaggtaggtggtcgcgcggttaacttgcttcctagtgaaaccgttcaggacatctacgggattgttgctaagaaagtcaacgagattctacaagcagacgcaatcaatgggaccgataacgaagtagttaccgtgaccgatgagaacactggtgaaatctctgagaaagtcaagctgggcactaaggcactggctggtcaatggctggcttacggtgttactcgcagtgtgactaagcgttcagtcatgacgctggcttacgggtccaaagagttcggcttccgtcaacaagtgctggaagataccattcagccagctattgattccggcaagggtctgatgttcactcagccgaatcaggctgctggatacatggctaagctgatttgggaatctgtgagcgtgacggtggtagctgcggttgaagcaatgaactggcttaagtctgctgctaagctgctggctgctgaggtcaaagataagaagactggagagattcttcgcaagcgttgcgctgtgcattgggtaactcctgatggtttccctgtgtggcaggaatacaagaagcctattcagacgcgcttgaacctgatgttcctcggtcagttccgcttacagcctaccattaacaccaacaaagatagcgagattgatgcacacaaacaggagtctggtatcgctcctaactttgtacacagccaagacggtagccaccttcgtaagactgtagtgtgggcacacgagaagtacggaatcgaatcttttgcactgattcacgactccttcggtaccattccggctgacgctgcgaacctgttcaaagcagtgcgcgaaactatggttgacacatatgagtcttgtgatgtactggctgatttctacgaccagttcgctgaccagttgcacgagtctcaattggacaaaatgccagcacttccggctaaaggtaacttgaacctccgtgacatcttagagtcggacttcgcgttcgcAtaa",
            gene_name="Pbad_rbs",
            first_six_aas="MDTINI"
        )
        self.rbsOptions.extend([opt1, opt2, opt3])

    def run(self, cds: str, ignores: set) -> RBSOption:
        """
        Selects an RBS that is not in the ignore set.
        
        Parameters:
            cds (str): The coding sequence.
            ignores (set): A set of RBS options to ignore.
        
        Returns:
            RBSOption: The selected RBS option.
        """
        for rbsopt in self.rbsOptions:
            if rbsopt not in ignores:
                return rbsopt
        raise Exception("No valid RBS option available.")

if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT..."

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)

'''