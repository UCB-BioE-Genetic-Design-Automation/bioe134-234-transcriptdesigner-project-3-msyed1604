class RNaseE_Checker:
    """
    Identifies potential RNase E cleavage sites within the designed RNA,
    which could lead to unintended RNA degradation.

    The function returns True if no cleavage sites are found.
    If a site is detected (False), it returns the problematic sequence.
    """

    def __init__(self):
        # Define common RNase E cleavage motifs (AU-rich sequences)
        self.cleavage_motifs = []

    def initiate(self):
        # Define common RNase E cleavage motifs (AU-rich sequences)
        self.cleavage_motifs = [
            "AUUUA", "UAUU", "AUUA", "UUAU", "AAUU", "UUAA"
        ]

    def run(self, rna_sequence: str) -> tuple[bool, str]:
        """
        Scans the RNA sequence for potential RNase E cleavage sites.
        
        Parameters:
            rna_sequence (str): The RNA sequence to check.
        
        Returns:
            tuple[bool, str]: A tuple where the first value is True if no cleavage 
                              sites are found, False otherwise. The second value 
                              is the problematic sequence if found.
        """
        # Convert to uppercase to avoid case sensitivity issues
        rna_sequence = rna_sequence.upper()

        # Scan for cleavage motifs
        for motif in self.cleavage_motifs:
            if motif in rna_sequence:
                return False, motif  # Cleavage site found

        return True, None  # No cleavage site found
