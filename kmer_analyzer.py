import sys #import python system module

def validate_sequence(sequence, k): #defines function validate_sequence that checks whether a DNA/RNA sequence is valid
    if len(sequence) < k: #if the length of the sequence is shorter than k, reject it
        return False
    for nucleotide in sequence:  
        if nucleotide in '1234567890':
            return False #loop through characters in the sequence, if any are a number, sequence is invalid, reject it
    return True #if sequence is k or longer and only letters, it passes, return true

def update_kmer_count(kmer_data, kmer, next_char): #makes function that is a dictionary that stores kmer statistics
    if kmer not in kmer_data: 
        kmer_data[kmer] = {'count': 1, 'next_chars': {}} #if a kmer has not been seen before, make an entry for it and store the number of times it appears as count, and the dictionary of what character follows after it 
    
    kmer_data[kmer]['count'] += 1  #increase count by one 
    
    if next_char not in kmer_data[kmer]['next_chars']: 
        kmer_data[kmer]['next_chars'][next_char] = 0 #if next character hasn't been see yet after kmer, initialize it
    kmer_data[kmer]['next_chars'][next_char] += 1 #increment how often the character follows kmer

    return kmer_data #return updated dictionary

def count_kmers_with_context(sequence, k): #builds full kmer data set for sequence
    kmer_data = {} #empty dictionary to store results
    
    for i in range(len(sequence) - k): #loop overall all possible kmer start positions, stop at the length of the sequence minus kmer length
        kmer = sequence[i:i+k] #extract out a substring of length k
        next_char = sequence[i+k] #this is the character that comes immediately after the kmer
        
        kmer_data = update_kmer_count(kmer_data, kmer, next_char) #update dictionary with observations
    
    return kmer_data #return kmer statistics


def write_results_to_file(kmer_data, output_filename): #define function to save kmer data to file
    sorted_kmers = sorted(kmer_data.keys()) #sorts kmers alphabetically 
    
    with open(output_filename, 'w') as f: #open output file for writing, overwrites existing file
        for kmer in sorted_kmers: #loops through each kmer 
            next_chars = kmer_data[kmer]['next_chars'] #get dictionary of next characters frequencies 
            
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            ) #builds formatted string, that writes out character and number of times observed, A:3, C:2 etc
            
            f.write(f"{kmer} {next_char_str}\n") #write one line per kmer 


def main(): 
    sequence_file = sys.argv[1] 
    k = int(sys.argv[2])
    output_file = sys.argv[3] #reads command line arguments of input file, kmer size, output file. So running in CL reads like python script.py DNAsequence.txt 4 output.file
    
    print(f"Reading sequences from {sequence_file}...") #print status message

    with open(sequence_file, 'r') as f: #opens input file, file with sequences
        for sequence in f:
            sequence = sequence.strip() #reads file line by line, which each line being a sequence, and remove newline character and spaces

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence") #if sequence is invalid, print warning and skip it
                continue
            
            kmer_data = count_kmers_with_context(sequence, k) #compute kmer statistics for sequence
            
            write_results_to_file(kmer_data, output_file) #write results to file

if __name__ == '__main__':
    main() #makes sure only runs when script is executed and not imported
