## About

Next-generation sequencing can currently produce hundreds of millions of reads
per lane of sample and that number increases at a dizzying rate.  Barcoding
individual sequences for multiple lines or multiple species is a cost-efficient
method to sequence and analyze a broad range of data.

Sabre is a tool that will demultiplex barcoded reads into separate files. 
It will work on both single-end and paired-end data in fastq format.
It simply compares the provided barcodes with each read and separates
the read into its appropriate barcode file, after stripping the barcode from
the read (and also stripping the quality values of the barcode bases).  If
a read does not have a recognized barcode, then it is put into the unknown file.  
Sabre also has an option (-m) to allow mismatches of the barcodes.

Sabre also supports gzipped file inputs.  Also, since sabre does not use the 
quality values in any way, it can be used on fasta data that is converted to
fastq by creating fake quality values.

Finally, after demultiplexing, sabre outputs a summary of how many records
went into each barcode file.
