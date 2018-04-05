## Requirements 

Sabre requires a C compiler; GCC or clang are recommended.  Sabre
relies on Heng Li's kseq.h, which is bundled with the source.

Sabre also requires Zlib, which can be obtained at
<http://www.zlib.net/>.


### Sabre Single End (`sabre se`)

`sabre se` takes an input fastq file and an input barcode data file and outputs 
the reads demultiplexed into separate files using the file names from the data file.
The barcodes will be stripped from the reads and the quality values of the barcode
bases will also be removed.  Any reads with unknown barcodes get put into the "unknown" 
file specified on the command line.  The -m option allows for mismatches in the barcodes.

#### Barcode data file format for single end

    barcode1 barcode1_output_file.fastq
    barcode2 barcode2_output_file.fastq
    etc...

Be aware that if you do not format the barcode data file correctly, sabre will not work properly.

#### Example

    sabre se -f input_file.fastq -b barcode_data.txt -u unknown_barcode.fastq
    sabre se -m 1 -f input_file.fastq -b barcode_data.txt -u unknown_barcode.fastq

### Sabre Paired End (`sabre pe`)

`sabre pe` takes two paired-end files and a barcode data file as input and outputs
the reads demultiplexed into separate paired-end files using the file names from the 
data file.  The barcodes will be stripped from the reads and the quality values of the barcode 
bases will also be removed.  Any reads with unknown barcodes get put into the "unknown" files 
specified on the command line.  It also has an option (-c) to remove barcodes from both files.  
Using this option means that if sabre finds a barcode in the first file, it assumes the paired 
read in the other file has the same barcode and will strip it (along with the quality values).  
The -m option allows for mismatches in the barcodes.

#### Barcode data file format for paired end

    barcode1 barcode1_output_file1.fastq barcode1_output_file2.fastq
    barcode2 barcode2_output_file1.fastq barcode2_output_file2.fastq
    etc...

Be aware that if you do not format the barcode data file correctly, sabre will not work properly.

#### Examples

    sabre pe -f input_file1.fastq -r input_file2.fastq -b barcode_data.txt \
    -u unknown_barcode1.fastq -w unknown_barcode1.fastq

    sabre pe -c -f input_file1.fastq -r input_file2.fastq -b barcode_data.txt \
    -u unknown_barcode1.fastq -w unknown_barcode1.fastq

    sabre pe -m 2 -f input_file1.fastq -r input_file2.fastq -b barcode_data.txt \
    -u unknown_barcode1.fastq -w unknown_barcode1.fastq
