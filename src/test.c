/* tab set at 4 spaces
 * using space characters, no tabs */

/* gcc -Wall -O3 -std=c99 -o test test.c -lz */

#include <stdio.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

int main (int argc, char *argv[]) {

    gzFile pe1=NULL;
    char *infn1=NULL;

    infn1 = (char*) malloc (strlen (argv[1]) + 1);
    strcpy (infn1, argv[1]);

    pe1 = gzopen (infn1, "r");
    kseq_t *fqrec1;
    fqrec1 = kseq_init (pe1);
    int l1;
    l1 = kseq_read (fqrec1);

    fprintf(stdout, "%d\n", l1);
    while ( (l1 = kseq_read (fqrec1) ) != -1) {
        //fprintf(stdout, "%s\n", fqrec1->name.s);
    }
    fprintf(stdout, "HERE\n");

    return 0;
    //gzFile pe2=NULL;
    //char *infn2=NULL;
    //infn2 = (char*) malloc (strlen (argv[2]) + 1);
    //strcpy (infn2, argv[2]);
    //pe2 = gzopen (infn2, "r");
    //kseq_t *fqrec2;
    //fqrec2 = kseq_init (pe2);
}
