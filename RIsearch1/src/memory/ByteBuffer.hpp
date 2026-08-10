#pragma once

#include <cstdlib>
#include <cstring>

/* A growable byte buffer that owns its storage.
 *
 * MallocRAII owns a block but cannot grow it, so callers that build a sequence
 * a character at a time end up managing the capacity themselves. This does that
 * part, and nothing else: no small-string optimisation, no allocator, no
 * exceptions, no templates.
 *
 * It grows with realloc, which lets the allocator extend the block in place
 * rather than allocate a second one and copy. That is the difference that
 * matters for a megabase record: std::string and std::vector both have to copy
 * at every doubling because their allocators cannot extend.
 */
class ByteBuffer {
public:
    ByteBuffer() = default;
    ~ByteBuffer() { std::free(m_data); }

    ByteBuffer(const ByteBuffer&) = delete;
    ByteBuffer& operator=(const ByteBuffer&) = delete;

    void clear() { m_size = 0; }

    [[nodiscard]] std::size_t size() const { return m_size; }
    [[nodiscard]] bool empty() const { return m_size == 0; }

    /* Valid only after terminate(). Empty rather than null before anything has
     * been written, so callers never have to check.
     */
    [[nodiscard]] const char* c_str() const { return m_data ? m_data : ""; }

    void append(const char* bytes, std::size_t count)
    {
        /* A blank line appends nothing, and memcpy may not be handed a null
         * destination even for zero bytes.
         */
        if (count == 0)
            return;
        reserve(m_size + count);
        std::memcpy(m_data + m_size, bytes, count);
        m_size += count;
    }

    void push_back(char byte)
    {
        reserve(m_size + 1);
        m_data[m_size++] = byte;
    }

    /* Writes the terminator without counting it, so c_str() can be handed to
     * anything expecting a C string while size() stays the byte count.
     */
    void terminate()
    {
        reserve(m_size + 1);
        m_data[m_size] = '\0';
    }

private:
    void reserve(std::size_t wanted)
    {
        if (wanted <= m_capacity)
            return;
        std::size_t capacity = m_capacity ? m_capacity : 128;
        while (capacity < wanted)
            capacity *= 2;
        m_data = static_cast<char*>(std::realloc(m_data, capacity));
        m_capacity = capacity;
    }

    char* m_data = nullptr;
    std::size_t m_size = 0;
    std::size_t m_capacity = 0;
};
