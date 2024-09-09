#this is a program to read in a fastq formatted file and exclude sequences
#that are shorter than or longer than a given length, as specified on line #38.
# the retained (non-excluded) sequences are written to the file "FILENAME_size_selected.fastq".
# the number of sequences in the new, size-sorted file will print on the command line

# Run the program on the command line with "...paa9$ python size_sort_fastq.py input_filename"

#!!! Output file should be "size_selected_fastq" check this.

## !!!! need to check Biopython usage, make sure it is installed (use Conda) also biopython.org!!!
## !!!! AND...check Nathan's Biopython slides, class #11.
import Bio.SeqIO
import sys



# get the sequence file name
seqfilename=sys.argv[1]

# get a base file name without all the extra and extensions
seqname = seqfilename.replace('_flt.fastq', '')

#open the FASTQ file and interate through all its sequences

seqfile= open(seqfilename)

counter=0
#!!!!need to check this line, not sure what the "w" specifies!!!! # SM: "w" stands for writing
out_handle=open(f"{seqname}_size_selected.fastq", "w")

for seq_record in Bio.SeqIO.parse(seqfile,"fastq"):

        length=len(seq_record.seq)
# trim last 4 bases off of sequence
        seq_record=seq_record[0:length-4]
# use new (trimmed) sequence length for sorting
        new_length=length-4
        if new_length >= 18 and new_length <= 24:
                counter=counter+1
                Bio.SeqIO.write(seq_record, out_handle,"fastq")

seqfile.close()
out_handle.close()

print(seqname)
print(counter)
