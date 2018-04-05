#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>
#include <time.h>
#include "sabre.h"
#include "utils.h"

    //more about getopts http://www.informit.com/articles/article.aspx?p=175771&seqNum=3
    static struct option paired_long_options[] = {
        {"fq1", required_argument, NULL, 'f'},
        {"fq2", required_argument, NULL, 'r'},
        {"barcodes", required_argument, NULL, 'b'},
        {"unassigned1", required_argument, NULL, 'u'},
        {"unassigned2", required_argument, NULL, 'w'},
        {"both-barcodes", optional_argument, NULL, 'c'},
        {"max-mismatch", required_argument, 0, 'm'},
        {"min-umi-len", required_argument, 0, 'l'},
        {"max-5prime-crop", required_argument, 0, 'a'},
        {"stats", required_argument, NULL, 's'},
        {"no-comment", no_argument, 0, 'n'},
        {"remove-seq", no_argument, 0, 'x'},
        {GETOPT_HELP_OPTION_DECL},
        {GETOPT_VERSION_OPTION_DECL},
        {NULL, 0, NULL, 0}
    };

void paired_usage (int status) {
    //TODO if --both-barcodes if given then --fq2 can't just be barcodes only 
    fprintf (stderr, "\n  Usage: %s pe [OPTIONS] -f <FASTQ_R1> -r <FASTQ_R2> -b <BARCODES>\
            \n\
            \n\
            \n  Options:\
            \n\
            \n    Required:\
            \n\
            \n        -f, --fq1 FILE            Input FASTQ R1 file, must have barcodes\
            \n        -r, --fq2 FILE            Input FASTQ R2 file\
            \n        -b, --barcodes FILE       Barcodes files, one barcode per line, e.g B\\tR1\\tR2\
            \n        -u, --unassigned1 FILE    Output unassigned R1 reads to [unassigned_R1.fastq]\
            \n        -w, --unassigned2 FILE    Output unassigned R2 reads to [unassigned_R1.fastq]\
            \n\
            \n    Other:\
            \n\
            \n        -c, --both-barcodes INT       Indicates that both FASTQ files have barcodes [0]\
            \n        -m, --max-mismatch INT        Maximum number of mismatches allowed in a barcode [0]\
            \n        -l, --min-umi-len INT         Minimum UMI length to keep [0]\
            \n        -a, --max-5prime-crop INT     Maximum number of possible bases cropped from 5prime [0]\
            \n        -n, --no-comment              Drop extra comments from FASTQ header [NULL]\
            \n        -x, --remove-seq              Replace --fq2 sequences with an N base. Only use if your R2 has just barcodes [NULL]\
            \n        -s, --stats FILE              Write stats to file instead of STDOUT [STDOUT]\
            \n\
            \n",
            PROGRAM_NAME);

    exit (status);
}

int paired_main (int argc, char *argv[]) {
    
    //clock_t begin = clock();
    time_t start, end;
    start = time(NULL);

    gzFile pe1=NULL;
    gzFile pe2=NULL;
    kseq_t *fqrec1;
    kseq_t *fqrec2;
    int l1,l2;
    FILE* barfile = NULL;
    FILE* unknownfile1=NULL;
    FILE* unknownfile2=NULL;
    FILE* log_file=NULL;
    int optc;
    extern char *optarg;
    char *infn1=NULL;
    char *infn2=NULL;
    char *barfn=NULL;
    char *unknownfn1=strdup("unassigned_R1.fastq");
    char *unknownfn2=strdup("unassigned_R2.fastq");
    int both_have_barcodes=0;
    barcode_data_paired *curr, *head, *temp;
    char barcode [MAX_BARCODE_LENGTH];
    char baroutfn1 [MAX_FILENAME_LENGTH];
    char baroutfn2 [MAX_FILENAME_LENGTH];
    int num_unknown=0;
    int total=0;
    int mismatch=0;

    int min_umi_len=0;
    int max_5prime_crop=0;
    char *log_fn=NULL;
    int no_comment=-1;
    int remove_seq=-1;

    while (1) {
        int option_index = 0;
        //colon after a flag means should have arguments and no colon means just a flag i.e bool, no args after it
        optc = getopt_long (argc, argv, "cnf:r:b:u:w:m:s:l:a:", paired_long_options, &option_index);

        if (optc == -1) break;

        switch (optc) {
            if (paired_long_options[option_index].flag != 0) break;

            case 'f':
            infn1 = (char*) malloc (strlen (optarg) + 1);
            strcpy (infn1, optarg);
            break;

            case 'r':
            infn2 = (char*) malloc (strlen (optarg) + 1);
            strcpy (infn2, optarg);
            break;

            case 'b':
            barfn = (char*) malloc (strlen (optarg) + 1);
            strcpy (barfn, optarg);
            break;

            case 'u':
            if(unknownfn1) {
                free(unknownfn1);
            }
            unknownfn1 = (char*) malloc (strlen (optarg) + 1);
            strcpy (unknownfn1, optarg);
            break;

            case 'w':
            if(unknownfn2) {
                free(unknownfn2);
            }
            unknownfn2 = (char*) malloc (strlen (optarg) + 1);
            strcpy (unknownfn2, optarg);
            break;

            case 'c':
            both_have_barcodes=1;
            break;

            case 'm':
            mismatch = atoi (optarg);
            break;

            case 's':
            log_fn = (char*) malloc (strlen (optarg) + 1);
            strcpy (log_fn, optarg);
            break;

            case 'l':
            min_umi_len = atoi (optarg);
            break;

            case 'a':
            max_5prime_crop = atoi (optarg);
            break;

            case 'n':
            no_comment = 1;
            break;

            case_GETOPT_HELP_CHAR(paired_usage);
            case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

            case '?':
            paired_usage (EXIT_FAILURE);
            break;

            default:
            paired_usage (EXIT_FAILURE);
            break;
        }
    }

    if (!infn1 || !infn2 || !unknownfn1 || !unknownfn2 || !barfn) {
        paired_usage (EXIT_FAILURE);
    }

    if (!strcmp (infn1, infn2) || !strcmp (infn1, unknownfn1) || !strcmp (infn1, unknownfn2) ||
            !strcmp (infn1, barfn) || !strcmp (infn2, unknownfn1) || !strcmp (infn2, unknownfn2) ||
            !strcmp (infn2, barfn) || !strcmp (unknownfn1, unknownfn2) || !strcmp (unknownfn1, barfn) ||
            !strcmp (unknownfn2, barfn)) {

        fprintf (stderr, "Error: Duplicate input and/or output file names.\n");
        return EXIT_FAILURE;
    }

    pe1 = gzopen (infn1, "r");
    if (!pe1) {
        fprintf (stderr, "Could not open input file 1 '%s'.\n", infn1);
        return EXIT_FAILURE;
    }

    pe2 = gzopen (infn2, "r");
    if (!pe2) {
        fprintf (stderr, "Could not open input file 2 '%s'.\n", infn2);
        return EXIT_FAILURE;
    }

    unknownfile1 = fopen (unknownfn1, "w");
    if (!unknownfile1) {
        fprintf (stderr, "Could not open unknown output file 1 '%s'.\n", unknownfn1);
        return EXIT_FAILURE;
    }

    unknownfile2 = fopen (unknownfn2, "w");
    if (!unknownfile2) {
        fprintf (stderr, "Could not open unknown output file 2 '%s'.\n", unknownfn2);
        return EXIT_FAILURE;
    }

    barfile = fopen (barfn, "r");
    if (!barfile) {
        fprintf (stderr, "Could not open barcode file '%s'.\n", barfn);
        return EXIT_FAILURE;
    }

    fprintf(stderr, "\n\
            \n  Running: %s\
            \n  Command line args:\
            \n      --fq1 %s\
            \n      --fq2 %s\
            \n      --barcodes %s\
            \n      --unassigned1 %s\
            \n      --unassigned2 %s\
            \n      --both-barcodes %d\
            \n      --max-mismatch %d\
            \n      --min-umi-len %d\
            \n      --max-5prime-crop %d\
            \n      --stats %s\
            \n      --no-comment %d\
            \n      --remove-seq %d\
            \n\
            \n  In Progess...\
            \n", PROGRAM_NAME,\
            infn1, infn2,\
            barfn,\
            unknownfn1, unknownfn2,\
            both_have_barcodes,\
            mismatch, min_umi_len, max_5prime_crop, log_fn, no_comment, remove_seq);

    /* Creating linked list of barcode data */
    // https://www.hackerearth.com/practice/data-structures/linked-list/singly-linked-list/tutorial/
    // where each node is represents one barcode from the barcode file
    // number of nodes should equal to number of barcodes (lines) in the file
    head = NULL;
    while (fscanf (barfile, "%s%s%s", barcode, baroutfn1, baroutfn2) != EOF) {
        curr = (barcode_data_paired*) malloc (sizeof (barcode_data_paired));
        curr->bc = (char*) malloc (strlen(barcode) + 1);
        strcpy (curr->bc, barcode);

        curr->bcfile1 = fopen (_mkdir(baroutfn1), "w");
        curr->bcfile2 = fopen (_mkdir(baroutfn2), "w");
        curr->num_records = 0;

        curr->next = head;
        head = curr;
    }

    fqrec1 = kseq_init (pe1);
    fqrec2 = kseq_init (pe2);
    // loop over all the reads and for every read loop over all barcodes and look for a match
    while ((l1 = kseq_read (fqrec1)) >= 0) {

        l2 = kseq_read (fqrec2);
        if (l2 < 0) {
            fprintf (stderr, "Error: PE file 2 is shorter than PE file 1. Disregarding rest of PE file 1.\n");
            break;
        }

        /* Go through all barcode data and check if any match to beginning of read */
        /* If it does then put read in that barcode's file, otherwise put in unknown file */
        curr = head;
        //NOTE here were are going through R1 only, that is -f,--fq1 reads and checking if they have a barcode
        //It appears that R1 is the file that must at least have barcodes and R2 may not have, but we are still going to 
        //appened it into FASTQ header of both R1 and R2 reads
        while (curr) {
            //zero means no mismatches found, that is barcode was found for that reads, therefore break and write it out
            if (strncmp_with_mismatch (curr->bc, fqrec1->seq.s, mismatch, max_5prime_crop) == 0) {
                break;
            }
            curr = curr->next;
        }

        if (curr != NULL) {
            // if UMI is shorter then 10, discard the reads
            if(strlen((fqrec1->seq.s)+strlen(barcode)) >= min_umi_len) {
                char *fqread1 = get_fqread(fqrec1, curr->bc, no_comment, -1);
                char *fqread2 = get_fqread(fqrec2, curr->bc, no_comment, remove_seq);

                fprintf(curr->bcfile1, "%s", fqread1);
                fprintf(curr->bcfile2, "%s", fqread2);

                free(fqread1);
                free(fqread2);

                curr->num_records += 2;
            }
        }

        else {
            char *fqread1 = get_fqread(fqrec1, NULL, no_comment, remove_seq);
            char *fqread2 = get_fqread(fqrec2, NULL, no_comment, remove_seq);
            
            fprintf(unknownfile1, "%s", fqread1);
            fprintf(unknownfile2, "%s", fqread2);

            free(fqread1);
            free(fqread2);

            num_unknown += 2;
        }

        total += 2;
    }

    if (l1 < 0) {
        l2 = kseq_read (fqrec2);
        if (l2 >= 0) {
            fprintf (stderr, "Error: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
        }
    }


    //if (!log_fn) { is this better?
    if (log_fn == NULL) {
        log_file = stdout;
    }
    else {
        log_file = fopen(log_fn, "w");
    }


    fprintf (log_file, "Barcode\tN_records\tN_pairs\tP_pairs\n");
    curr = head;
    int total_pairs = total/2; 

    while (curr) {

        int n_pairs = curr->num_records/2;
        float percent_pairs = (float) n_pairs/total_pairs;

        fprintf (log_file,"%s\t%d\t%d\t%.2f\n", curr->bc, curr->num_records, n_pairs, percent_pairs);
        curr = curr->next;
    }

    int unknown_pairs = num_unknown/2;
    float percent_unknown = (float) unknown_pairs/total_pairs;
    float tot_chk = (float) total_pairs/total_pairs;

    fprintf (log_file, "unassigned\t%d\t%d\t%.2f\n", num_unknown, unknown_pairs, percent_unknown);
    fprintf (log_file, "total\t%d\t%d\t%.2f\n", total, total_pairs, tot_chk);

    end = time(NULL);
    fprintf(stderr, "\n All done :) \
	             \n It took %.2f minutes\n", difftime(end, start)/60);

    kseq_destroy (fqrec1);
    kseq_destroy (fqrec2);
    gzclose (pe1);
    gzclose (pe2);
    fclose (unknownfile1);
    fclose (unknownfile2);
    fclose (barfile);
    fclose (log_file);

    free (infn1);
    free (infn2);
    free (barfn);
    free (unknownfn1);
    free (unknownfn2);
    free (log_fn);

    curr = head;
    while (curr) {
        fclose (curr->bcfile1);
        fclose (curr->bcfile2);
        free (curr->bc);
        temp = curr;
        curr = curr->next;
        free (temp);
    }

    return EXIT_SUCCESS;
}
