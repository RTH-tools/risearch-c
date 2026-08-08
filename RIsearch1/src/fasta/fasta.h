/* fasta.h
 * Declarations for simple FASTA i/o library
 * SRE, Sun Sep  8 05:37:38 2002 [AA2721, transatlantic]
 * CVS $Id$
 */

#pragma once

#include <cstdio>
#include <cstdint>

#define FASTA_MAXLINE 512 /* Requires FASTA file lines to be <512 characters */

typedef struct fastafile_s {
    FILE* fp;
    char buffer[FASTA_MAXLINE];
} FASTAFILE;

FASTAFILE* OpenFASTA(const char* seqfile);
int ReadFASTA(FASTAFILE* fp, char** ret_seq, char** ret_name, std::uint32_t* ret_L);
void CloseFASTA(FASTAFILE* ffp);

