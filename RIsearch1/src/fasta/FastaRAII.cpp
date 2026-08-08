#include "FastaRAII.h"
#include "fasta.h"


FastaRAII::FastaRAII(const char* sequence_file) : m_sequence_file(sequence_file)
{
    m_fasta_file_handle = OpenFASTA(m_sequence_file);
}

FASTAFILE* FastaRAII::handle() const
{
    return m_fasta_file_handle;
}

FastaRAII::~FastaRAII()
{
    if (m_fasta_file_handle) {
        CloseFASTA(m_fasta_file_handle);
    }
}
