import re

# --- PART 1: GENETIC MAPPING ---
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
    'TGA':'_', 'GTC':'V', 'GTT':'V', 'GTA':'V', 'GTG':'V', 'TGC':'C', 'TGT':'C', 'TGG':'W',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q'
}

# --- PART 2: PROCESSING ENGINE ---

def clean_sequence(raw_data: str) -> str:
    """Filters non-genomic noise and standardizes sequence to uppercase using regex."""
    clean_data = raw_data.upper()
    return re.sub(r'[^ATGC]', '', clean_data)

def translate_dna(dna_seq: str, stop_at_stop_codon: bool = False) -> str:
    """Translates DNA triplets into primary protein amino acid structure efficiently."""
    protein = []
    for i in range(0, (len(dna_seq) // 3) * 3, 3):
        codon = dna_seq[i:i+3]
        aa = GENETIC_CODE.get(codon, "?")
        if stop_at_stop_codon and aa == '_':
            break
        protein.append(aa)
    return "".join(protein)

# --- PART 3: MOLECULAR DIAGNOSTIC ENGINE ---

def run_hbb_diagnostic(protein_seq: str) -> str:
    """Specific biomarker detection for Sickle Cell Anemia (HbS)."""
    if len(protein_seq) >= 6:
        target_aa = protein_seq[5]  # Index 5 is Position 6
        report = ["\n" + "*"*50, " [MOLECULAR DIAGNOSTIC REPORT: HBB GENE] ", "*"*50]
        
        if target_aa == 'V':
            report.append("FINAL STATUS: POSITIVE (Pathogenic)")
            report.append("FINDINGS: HbS mutation (Valine) detected at Position 6.")
        elif target_aa == 'E':
            report.append("FINAL STATUS: NEGATIVE (Normal)")
            report.append("FINDINGS: Glutamic Acid (Wild Type) detected at Position 6.")
        else:
            report.append(f"FINAL STATUS: ATYPICAL VARIANT ('{target_aa}')")
            report.append("WARNING: Non-standard amino acid detected. Clinical review required.")
        report.append("*"*50)
        
        result = "\n".join(report)
        print(result)
        return result
    else:
        err = "\nCRITICAL ERROR: DNA sequence length is insufficient for HBB analysis."
        print(err)
        return err

# --- PART 4: INTERACTIVE INTERFACE ---

def main_console():
    print("\n" + "="*55)
    print("      GENOMIC ANALYSIS SUITE - V3.1 PROFESSIONAL")
    print("="*55)
    print("1. PROTEIN TRANSLATION (General)")
    print("2. SICKLE CELL ANEMIA DETECTION (HBB Gene)")
    print("3. SEQUENCE METRICS (GC Content & Length)")
    print("4. TERMINATE PROTOCOL")
    
    user_choice = input("\nSelect system option (1-4): ")

    if user_choice == '1':
        print("\n--- GENERAL TRANSLATION ENGINE ---")
        user_dna = input("Enter DNA sequence (Leave empty for Human Insulin sample): ").strip()
        if not user_dna:
            user_dna = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCCCTGCTGGCCCTCTGGGGACCTGACCCAGCC"
            print("Using default Human Insulin sequence...")
        
        clean_dna = clean_sequence(user_dna)
        print(f"\nPRIMARY STRUCTURE OUTPUT: {translate_dna(clean_dna)}")
        
    elif user_choice == '2':
        print("\n--- INITIATING HBB SCANNER ---")
        print("Input DNA sequence for HBB analysis:")
        dna_input = input("DNA SEQUENCE > ")
        clean_dna = clean_sequence(dna_input)
        prot_output = translate_dna(clean_dna)
        print(f"\nTRANSLATION RESULT: {prot_output}")
        run_hbb_diagnostic(prot_output)
        
    elif user_choice == '3':
        dna_input = input("\nEnter sequence for metric calculation: ")
        clean_dna = clean_sequence(dna_input)
        if clean_dna:
            gc_val = (clean_dna.count("G") + clean_dna.count("C")) / len(clean_dna) * 100
            print("\n--- DATA METRICS ---")
            print(f"GC CONTENT: {gc_val:.2f}%")
            print(f"TOTAL BASES: {len(clean_dna)} bp")
        else:
            print("\nERROR: No valid genomic input detected.")

    elif user_choice == '4':
        print("\nSystem offline. Session ended.")
        return

    main_console()

if __name__ == "__main__":
    main_console()
