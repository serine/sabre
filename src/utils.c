#include "utils.h"

// https://stackoverflow.com/questions/2336242/recursive-mkdir-system-call-on-unix/11425692
// https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c
const char * _mkdir(const char *file_path) {

    if(!file_path) {
        fprintf (stderr,
                "ERROR: This shouldn't happend, file path == %s\n",
                file_path);
        exit(EXIT_FAILURE);
    }

    // return straigth away if a file_path is not nested file path
    if(strstr(file_path, "/") == NULL) {
        return file_path;
    }
    //TODO check if directory already exists or not
    //struct stat st = {0};
    //
    //if (stat(dir, &st) == -1) {
    //    mkdir(tmp, S_IRWXU);
    //}

    //char tmp[PATH_MAX]; // can't get this to work..
    char tmp[256];
    char *p = NULL;
    char *dirc = strdup(file_path);
    size_t len;

    const char *dir = dirname(dirc);
    snprintf(tmp, sizeof(tmp),"%s", dir);
    len = strlen(tmp);

    if(tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }

    for(p = tmp + 1; *p; p++) {
        if(*p == '/') {
            *p = 0;
            mkdir(tmp, S_IRWXU);
            *p = '/';
        }
    }

    mkdir(tmp, S_IRWXU);
    free(dirc);

    return file_path;
}

//NOTE retuns zero on success
//strcmp can be used for sorting, returns pos, zero, neg
//BUT this new implementation can't be used as such just FYI
int chk_bc_mtch(const char *orig_bc, const char *orig_read, size_t mismatch, int max_5prime_crop) {
    int orig_read_len = strlen(orig_read);
    int orig_bc_len = strlen(orig_bc);
    int n_crop = 0;
    if(orig_bc_len > orig_read_len) {
        fprintf (stderr,
                "WARNING: Length of the barcode %d is greater than length of the reads %d.",
                orig_bc_len, orig_read_len);
        return -1;
    }

    while(n_crop <= max_5prime_crop) {

        if(n_crop > orig_read_len) {
            return -1;
        }

        int cnt = 0;
        char u1, u2;
        const char *bc = orig_bc;
        const char *read = orig_read+n_crop;
        int bc_len = orig_bc_len;

        while (bc_len-- > 0) {
            u1 = *bc++;
            u2 = *read++;

            if (u1 != u2) {
                cnt++;
                if (cnt > mismatch) {
                    break;
                }
            }

            if (u1 == '\0' || u2 == '\0') {
                return n_crop;
            }
        }

        if(cnt <= mismatch) {
            return n_crop;
        }

        n_crop++;
    }
    //this is in the case of error
    return -1;
}

// https://stackoverflow.com/questions/21880730/c-what-is-the-best-and-fastest-way-to-concatenate-strings
//TODO this is a fastq mystrcat function, that returns a pointer to the end of the string
void get_fqread(char *fqread, fq_rec_t *fq_rec, char *barcode, char *umi_idx, int no_comment, int n_crop) {

    fqread[0] = '\0';

    if(n_crop < 0) {
        fprintf(stderr,
                "ERROR: n_crop set to %d. This can't happend\n",
                n_crop);
        exit(EXIT_FAILURE);
    }

    //@READNAME:BACRCODE:UMI
    //1st line
    strcat(fqread, fq_rec->name);
    //TODO later can have conditional here depending on the the structure and/or BARCODE/UMI
    if(barcode) {
        strcat(fqread, ":");
        strcat(fqread, barcode);
    }
    else if(!barcode) {
        barcode = "";
    }

    if(umi_idx) {
        strcat(fqread, ":");
        strcat(fqread, umi_idx);
    }

    if(fq_rec->comment && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fq_rec->comment);
    }
    strcat(fqread, "\n");

    //2nd line
    strcat(fqread, (fq_rec->seq)+strlen(barcode)+n_crop);
    strcat(fqread, "\n");

    //3rd line
    strcat(fqread, "+");
    strcat(fqread, fq_rec->name);
    if(fq_rec->comment && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fq_rec->comment);
    }
    strcat(fqread, "\n");

    //4th line
    strcat(fqread, (fq_rec->qual)+strlen(barcode)+n_crop);
    strcat(fqread, "\n");
}

void get_merged_fqread(char *fqread, fq_rec_t *fq_rec1, fq_rec_t *fq_rec2, char *barcode, char *umi_idx, int no_comment, int n_crop) {
    fqread[0] = '\0';
    //@READNAME:BACRCODE:UMI
    //1st line
    strcat(fqread, fq_rec1->name);
    //TODO later can have conditional here depending on the the structure and/or BARCODE/UMI
    if(barcode) {
        strcat(fqread, ":");
        strcat(fqread, barcode);
    }

    if(umi_idx) {
        strcat(fqread, ":");
        strcat(fqread, umi_idx);
    }

    if(fq_rec1->comment && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fq_rec1->comment);
    }
    strcat(fqread, "\n");

    //2nd line
    strcat(fqread, fq_rec2->seq);
    strcat(fqread, "\n");

    //3rd line
    strcat(fqread, "+");
    strcat(fqread, "");
    //strcat(fqread, fq_rec2->name);
    //if(fq_rec2->comment && no_comment == -1) {
    //    strcat(fqread, " ");
    //    strcat(fqread, fq_rec2->comment);
    //}
    strcat(fqread, "\n");

    //4th line
    strcat(fqread, (fq_rec2->qual));
    strcat(fqread, "\n");
}

void get_bc_fn(char **bcout_fn, char *s_name, char *barcode, int read_type) {

    if(strlen(s_name) > MAX_FILENAME_LENGTH) {
        fprintf (stderr,
                "ERROR: Too many characters in your sample name; %s:%zd \n",
                s_name, strlen(s_name));
        exit(EXIT_FAILURE);
    }

    strcat(*bcout_fn, s_name);
    strcat(*bcout_fn, "_");
    strcat(*bcout_fn, barcode);

    if(read_type == 1) {
        strcat(*bcout_fn, "_R1.fastq");
    }
    else if(read_type == 2) {
        strcat(*bcout_fn, "_R2.fastq");
    }
    else {
        fprintf (stderr,
                "ERROR: This shouldn't happen, wrong read type was passed through -> %d\n",
                read_type);
        exit(EXIT_FAILURE);
    }
}

void do_rev_comp(char *seq) {
    // this code was taken from seqtk with slight modification to drop kseq struct
    // https://github.com/lh3/seqtk
    char comp_tab[] = {
         0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
        16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
        32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
        48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
    };

    int seq_len = strlen(seq);
    int c0, c1;
    // reverse complement sequence
    for (int i = 0; i < seq_len>>1; ++i) {
        c0 = comp_tab[(int)seq[i]];
        c1 = comp_tab[(int)seq[seq_len - 1 - i]];
        seq[i] = c1;
        seq[seq_len - 1 - i] = c0;
    }
    // complement the remaining base
    if(seq_len & 1) {
        seq[seq_len>>1] = comp_tab[(int)seq[seq_len>>1]];
    }
    //if (fqrec1->qual.l) {
    //    for (i = 0; i < fqrec1->seq.l>>1; ++i) // reverse quality
    //        c0 = fqrec1->qual.s[i], fqrec1->qual.s[i] = fqrec1->qual.s[fqrec1->qual.l - 1 - i], fqrec1->qual.s[fqrec1->qual.l - 1 - i] = c0;
    //}
}

void set_default_params(param_t *params) {
    params->mismatch = 0;
    params->combine = -1;
    params->umi = -1;
    params->paired = -1;
    params->min_umi_len = 0;
    params->max_5prime_crop = 0;
    params->no_comment = -1;
}

void params_destroy(param_t *params) {
    gzclose(params->fq1_fd);
    gzclose(params->fq2_fd);
    fclose(params->unassigned1_fd);
    fclose(params->unassigned2_fd);
    fclose(params->umis_2_short_fd);
}
