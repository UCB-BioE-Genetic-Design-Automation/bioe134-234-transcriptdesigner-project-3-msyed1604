import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence using a hybrid 
    Guided Random + Sliding Window approach for optimal codon selection satisfying
    high CAI, low hairpin count, and the absence of internal promoters & forbidden sequences.
    """


    def initiate(self):
        """
        Initialization method for RBSChooser
        """

        self.aminoAcidToCodon = {}
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.promoterChecker = PromoterChecker()
        self.codonChecker = CodonChecker()
        self.rbsChooser = RBSChooser()
        

         # Initialize checkers
        self.forbiddenChecker.initiate()
        self.promoterChecker.initiate()
        self.codonChecker.initiate()
        self.rbsChooser.initiate()

        # Load codon usage data from the provided text file
        self.aminoAcidToCodon = self.load_codon_usage('genedesign/data/codon_usage.txt')

    def load_codon_usage(self, filepath: str) -> dict:
        """
        Parses a codon usage file and returns a dictionary mapping amino acids to their codons and frequencies.
        
        Parameters:
            filepath (str): Path to the codon usage file.
        
        Returns:
            dict: A dictionary where keys are amino acids and values are lists of tuples (codon, frequency).
        """
        amino_acid_to_codon = {}

        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
            
                if len(parts) >= 3:
                    codon = parts[0].strip()  # Codon (e.g., TTT)
                    aa = parts[1].strip()     # Amino acid (e.g., F)
                    frequency = float(parts[2])  # Frequency (e.g., 0.58)

                # Add codon and frequency to the corresponding amino acid entry
                    if aa not in amino_acid_to_codon:
                        amino_acid_to_codon[aa] = []
                    amino_acid_to_codon[aa].append((codon, frequency))
        
        return amino_acid_to_codon
    
    def guided_random(self, peptide: str) -> list:
        """
        Generates an initial random DNA sequence using codon frequencies from the loaded data.
        
        Parameters:
            peptide (str): The protein sequence to translate.
        
        Returns:
            List[str]: A list of selected codons representing the DNA sequence.
        """
        dna_sequence = []

        for aa in peptide:
            # Get possible codons and their frequencies for this amino acid
            possible_codons = self.aminoAcidToCodon.get(aa)
            if possible_codons:
                # Extract codons and their corresponding weights (frequencies)
                codons, weights = zip(*possible_codons)
                # Choose a random codon based on its frequency
                selected_codon = random.choices(codons, weights=weights)[0]
                dna_sequence.append(selected_codon)

        return dna_sequence

    def validate_window(self, window: list) -> bool:
        """
        Validates a 3-codon window by checking it against all validation criteria.
        
        Parameters:
            window (list): A 3-codon window to validate.
        
        Returns:
            bool: True if the window passes all checks; False otherwise.
        """
        
        # Convert window into DNA sequence string
        dna_seq = ''.join(window)

        # Run all checkers on this window
        
        # Forbidden Sequence Check
        forbidden_okay = self.forbiddenChecker.run(dna_seq)
        
        # Hairpin Check
        hairpin_okay = hairpin_checker(dna_seq)
        
        # Promoter Check
        promoter_okay = self.promoterChecker.run(dna_seq)
        
        # Codon Quality Check (diversity, rare count, CAI)
        quality_okay = self.codonChecker.run(window)

        return forbidden_okay and hairpin_okay and promoter_okay and quality_okay
    
    def sliding_window_refinement(self, codons: list) -> list:
        """
        Refines a DNA sequence using a sliding window approach with validation using checkers.
        
        Parameters:
            codons (list): The initial list of codons to refine.
        
        Returns:
            List[str]: The refined list of codons after applying sliding window optimization and validation.
        """

        # This will act as the preamble (upstream sequence)
        refined_codons = [] 
         
        # Iterate over windows of 3 amino acids (9 nucleotides)
        for i in range(0, len(codons) - 2):
            # Current 3-codon window
            window = codons[i:i+3]

            # Optimize this window based on upstream and downstream context
            optimized_window = self.optimize_window(window)
            
            # Keep only middle codon from optimized window
            refined_codons.append(optimized_window[1])
        
        # Append last two remaining codons after final window
        refined_codons.extend(codons[-2:])
        
        return refined_codons

    def optimize_window(self, window: list) -> list:
        """
        Optimizes a 3-codon window by trying alternative middle codons until validation passes.
        
        Parameters:
            window (list): A 3-codon window to optimize.
        
        Returns:
            List[str]: The optimized 3-codon window that passes validation.
        
         Raises:
             Exception: If no valid alternative is found after multiple attempts.
         """
        
        max_attempts = 10  # Limit attempts to find valid alternatives
        
        for _ in range(max_attempts):
            if self.validate_window(window):
                return window  # If valid return this window
            
             # Try another random alternative for this window's middle codon
            middle_aa = [self.aminoAcidToCodon[window[1][0]]]
            possible_codons = [random.choice(middle_aa)]
            
            if possible_codons:
                new_middle_codon = random.choice(possible_codons)[0]
                window[1] = new_middle_codon

        raise Exception("Failed to find valid alternative after multiple attempts")

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA using hybrid Guided Random + Sliding Window approach,
         and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        
        # Step 1: Generate initial random design
        initial_codons = self.guided_random(peptide)
        
        # Step 2: Refine using sliding window optimization
        refined_codons = self.sliding_window_refinement(initial_codons)

        # Step 3: Append stop codon
        refined_codons.append("TAA")

        # Step 4: Build CDS from refined codons
        cds = ''.join(refined_codons)

        # Step 5: Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Step 6: Return Transcript object
        return Transcript(selectedRBS, peptide, refined_codons)


if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    
    # Provide path to your 'codon_usage.txt' file here
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    print(transcript)

'''
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # Translate peptide to codons
        codons = [self.aminoAcidToCodon[aa] for aa in peptide]

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)

'''