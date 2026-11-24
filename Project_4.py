############ Project 4 ############ 

import gzip
import math
import argparse
import matplotlib.pyplot as plt

converge_threshold=0.00001 # |H_current-H_previous|
check_interval=1000 #compute the Shannon entropy every N reads

def read_fastq(filename): #reads a fastq.gz and gives the sequences one by one
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i%4==1: #we get only the lines with sequences
                yield line.strip()


def find_kmers(sequence, k): #gets all kmers from a given sequence
    sequence=sequence.upper()
    for i in range(len(sequence)-k+1): #last possible position for a kmer to start
        kmer=sequence[i:i+k]
        if 'N' not in kmer:
            yield kmer


def calculate_shannon_entropy(kmer_counts): #calculates entropy from a dictionary with the kmers counts
    total_kmers=sum(kmer_counts.values())
    if total_kmers==0:
        return 0.0
    entropy=0.0
    for count in kmer_counts.values():
        p_i=count/total_kmers
        entropy-=p_i*math.log2(p_i)
    return entropy
#Low entropy means the sample has almost one kmer
#High entropy means the sample haw many different kmers and all of them have similar frequency 


def plot_richness(processed_reads, unique_kmers, k, output_name): #processed_reads:list with the check intervals, unique_kmers:list with number of kmers in each interval
    plt.figure(figsize=(10,6))
    plt.plot(processed_reads, unique_kmers, marker='.', linestyle='-')
    plt.title(f"Metagenome Richness: Unique k-mers (k={k}) pre processed reads", fontsize=14)
    plt.xlabel("Reads Processed", fontsize=12)
    plt.ylabel("Unique k-mers Discovered", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_name)
    print(f"\n(k={k}) Plot saved in: {output_name}")
#climbing line at the start: we find new kmers (new species)
#horizontal line: we stop finding new kmers, we have seen the most diversity, the entropy has converged
#climbing line all the way: very rich sample, we have not seen the most diversity even if the entropy has converged


def process_file(filename, k):
    kmer_counts={} #dictionary to count kmers
    reads_processed=0 #how many DNA sequences we read
    total_kmers_found=0 
    current_entropy=0.0 #the entropy we calculated on the last interval
    previous_entropy=0.0 #the entropy on the previous interval
    reads_list=[] #list with reads for the plot
    unique_kmers_list=[] #list with unique kmers for the plot

    sequence_generator=read_fastq(filename) #gives sequences one by one
    convergence_reached=False #we stop when we have convergence

    for sequence in sequence_generator: #for each DNA sequence (we get from the generator() of the file
        reads_processed+=1 

        for kmer in find_kmers(sequence, k): #we get all kmers in the sequence
            kmer_counts[kmer]=kmer_counts.get(kmer,0)+1 #we add the kmer in the dictionary and keep the counts of it
            total_kmers_found+=1

        #we check for convergence in every interval
        if reads_processed%check_interval==0:
            previous_entropy=current_entropy
            current_entropy=calculate_shannon_entropy(kmer_counts) #we calculate the entropy with the kmers we have until now
            unique_kmers=len(kmer_counts) #how many unique kmers we have

            reads_list.append(reads_processed) #x-axis of the plot
            unique_kmers_list.append(unique_kmers) #y-axis of the plot

            difference=abs(current_entropy-previous_entropy)

            if difference<converge_threshold and reads_processed>check_interval: #we want to reach at least the second interval
                convergence_reached=True
                break
    
    #after the for loop is finished (either we reached convergence or there are no more sequences in the file)
    if not convergence_reached:
        print(f"\n(k={k}) The file is finished after {reads_processed} reads")
        if reads_processed%check_interval!=0: #we calculate the entropy after the last read
            current_entropy=calculate_shannon_entropy(kmer_counts)
            unique_kmers=len(kmer_counts)
            reads_list.append(reads_processed)
            unique_kmers_list.append(unique_kmers)

    #either we have convergence or not
    unique_kmers=len(kmer_counts)
    plot_filename=f"Richness_plot_for_k={k}.png"

    #we make the plot
    plot_richness(reads_list, unique_kmers_list, k, plot_filename)

    return {'k':k, 'shannon_entropy':current_entropy, 'reads_processed':reads_processed, 'total_kmers':total_kmers_found, 'unique_kmers':unique_kmers, 'plot_file':plot_filename} #we return a dictionary with all results
    


parser=argparse.ArgumentParser()

parser.add_argument("fastq_file")
parser.add_argument("k", type=int, nargs='?', default=None) #? means it is optional

args=parser.parse_args()

all_results=[] #list for all results

if args.k is not None:
    if not (4<=args.k<=7):
        print("Error: k must be between 4-7")
        exit(1) #it stops the script
    
    result=process_file(args.fastq_file, args.k)
    if result:
        all_results.append(result)
else: #if k=None (Multi-k)
    for k in range(4,8):
        result=process_file(args.fastq_file, k)
        if result:
            all_results.append(result)

if all_results:
    for result in all_results:
        print(f"k={result['k']},    Shannon Entropy={result['shannon_entropy']},    Reads Processed={result['reads_processed']},    Total k-mers={result['total_kmers']},   Unique k-mers={result['unique_kmers']},     Plot file={result['plot_file']}")
else:
    print("Analysis failed")