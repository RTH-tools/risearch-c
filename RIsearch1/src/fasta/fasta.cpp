/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id$
 */


/* Function: OpenFASTA(), ReadFASTA(), CloseFASTA().
 * Date:     SRE, Sun Sep  8 06:39:26 2002 [AA2721, transatlantic]
 *
 * Purpose:  A very rudimentary FASTA file reading API. Designed
 *           for simplicity and clarity, not for robustness.
 *
 *           The API is:
 *
 *           ffp = OpenFASTA(seqfile);
 *           while (ReadFASTA(ffp, seq, name))
 *           {
 *             do stuff with sequence;
 *           }
 *           CloseFASTA(ffp);
 *
 * Args:
 *           seqfile   - name of a FASTA file to open.
 *           seq       - RETURN: one sequence
 *           name      - RETURN: name of the sequence
 *           ffp       - ptr to a FASTAFILE object.
 *
 * Commentary:
 *           The basic problem with reading FASTA files is that there is
 *           no end-of-record indicator. When you're reading sequence n,
 *           you don't know you're done until you've read the header line
 *           for sequence n+1, which you won't parse 'til later (when
 *           you're reading in the sequence n+1). One common trick for
 *           this is to implement a one-line "lookahead" buffer that you
 *           can peek at, before parsing later.
 *
 *           This buffer is kept in a small structure (a FASTAFILE), rather
 *           than in a static char[] in the function. This allows
 *           us to have multiple FASTA files open at once. The static approach
 *           would only allow us to have one file open at a time. ANSI C
 *           predates the widespread use of parallel programming. It was
 *           not overly concerned about the drawbacks of statics. Today,
 *           though, you should keep in mind that you may someday want to
 *           turn your program into a multithreaded, parallel program, and
 *           all functions in parallelized code must be "reentrant": able to
 *           be called a second time - with different arguments,
 *           and while the code in the first function call is still executing! -
 *           without overwriting or corrupting any static storage in the
 *           function. Statics have fewer uses now (for example, to
 *           test that some initialization code for a function is run once
 *           and only once.)
 *
 * Limitations:
 *           There is no error handling, for clarity's sake. Also,
 *           the parser is brittle. Improper FASTA files (for instance,
 *           blank lines between records) will cause unexpected
 *           behavior. Real file parsers are more complex.
 *           In real life, they have to deal with absolutely anything the user might
 *           pass as a "FASTA file"; and either parse it correctly,
 *           or detect that it's an invalid format and fail cleanly.
 *
 *           Lines are read in from the file using ANSI C's fgets(). fgets()
 *           requires a maximum buffer length (here, FASTA_MAXLINE, which is
 *           defined as 512 in bio5495.h). Some FASTA files have very long
 *           description lines, however; notably the NCBI NR database. Static
 *           limitations on things like line or sequence lengths should be
 *           avoided. An example of a replacement for fgets() that dynamically
 *           allocates its buffer size and allows any line length is
 *           SQUID's sre_fgets().
 *
 *           We use ANSI C's strtok() to parse the sequence name out of the line.
 *           strtok() is deprecated in modern programs because it is not threadsafe.
 *           (See comments above.) An example of a threadsafe version is
 *           SQUID's sre_strtok().
 *
 * Returns:
 *           OpenFASTA() returns a FASTAFILE pointer, or NULL on failure (for
 *           instance, if the file doesn't exist, or isn't readable).
 *
 *           ReadFASTA() returns true on success, or false if there are no
 *           more sequences to read in the file.
 *
 *           CloseFASTA() "always succeeds" and returns void.
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "fasta.h"
#include "fasta/ResidueTable.h"


namespace {

/* Besides residues, a sequence line may hold the line ending, the spaces and
 * digits of a coordinate column, and the '-' and '.' of an alignment gap. A
 * byte outside that set is corruption, and dropping it quietly would return a
 * shorter sequence that still looks valid.
 */
bool tolerated_in_sequence(unsigned char c)
{
    return c == '\n' || c == '\r' || c == ' ' || c == '\t' || c == '-' || c == '.' ||
           (c >= '0' && c <= '9');
}

/* Out of line so that fprintf does not sit inside the per-character loop, where
 * it costs about a tenth of the throughput for a branch a readable file never
 * takes.
 */
[[noreturn]] __attribute__((noinline, cold)) void reject_corrupt_byte(unsigned char c,
                                                                     const char* name)
{
    fprintf(stderr, "Corrupt byte 0x%02x in the sequence of '%s'.\n", c, name);
    exit(1);
}

} // namespace


FASTAFILE* OpenFASTA(const char* seqfile)
{
    FASTAFILE* ffp = reinterpret_cast<FASTAFILE*>(malloc(sizeof(FASTAFILE)));

    if (strcmp(seqfile, "-")) {        /*returns 0/FALSE if they are same! */
        ffp->fp = fopen(seqfile, "r"); /* Assume seqfile exists & readable!   */
    } else {
        ffp->fp = stdin;
    }
    if (ffp->fp == NULL) {
        free(ffp);
        return NULL;
    }
    if ((fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) == NULL) {
        free(ffp);
        return NULL;
    }
    return ffp;
}

bool ReadFASTA(FASTAFILE* ffp, ByteBuffer& ret_seq, ByteBuffer& ret_name)
{
    /* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline.
     */
    if (ffp->buffer[0] != '>')
        return false;

    /* Parse out the name: the first non-whitespace token after the >
     */
    const char* s = strtok(ffp->buffer + 1, " \t\r\n");
    ret_name.clear();
    ret_name.append(s, strlen(s));
    ret_name.terminate();

    /* Everything else 'til the next descline is the sequence. clear() keeps the
     * capacity the previous record grew, so reading a file of similar records
     * settles on one buffer instead of allocating per record.
     */
    ret_seq.clear();
    while (fgets(ffp->buffer, FASTA_MAXLINE, ffp->fp)) {
        if (ffp->buffer[0] == '>')
            break; /* a-ha, we've reached the next descline */

        /* A sequence line is normally residues followed by a line ending, so
         * take the leading run in one append and only go character by character
         * over whatever follows it.
         */
        const unsigned char* line = reinterpret_cast<const unsigned char*>(ffp->buffer);
        std::size_t run = 0;
        while (kResidue.is[line[run]])
            run++;
        ret_seq.append(ffp->buffer, run);

        for (std::size_t i = run; line[i] != '\0'; i++) {
            const unsigned char c = line[i];
            if (!kResidue.is[c]) {
                if (!tolerated_in_sequence(c))
                    reject_corrupt_byte(c, ret_name.c_str());
                continue; /* accept any alphabetic character */
            }
            ret_seq.push_back(static_cast<char>(c));
        }
    }
    ret_seq.terminate();
    return true;
}

void CloseFASTA(FASTAFILE* ffp)
{
    fclose(ffp->fp);
    free(ffp);
}


/* what follows is a useful idiom: when you're writing a .c file that's supposed
 * to be a module of library functions, include one or more "test drivers".
 * These are small main()'s, normally ifdef'ed out of the code, that
 * enable the .c file to be compiled into one or more standalone test programs.
 * This lets you test your module in relative isolation, which tends
 * to lead to faster debugging and more robust code. It also
 * provides a convenient way to document a working minimal API: for example,
 * the main() here is a minimal FASTA reader. And it also tends to
 * have a useful psychological effect on you: it tends to encourage you
 * to simplify your APIs, so that small test programs can demonstrate
 * the full power of the API.
 */
