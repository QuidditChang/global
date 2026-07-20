/* Minimal NumPy-compatible, uncompressed NPZ writer. */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "npz_writer.h"

static void put_u16(FILE *file, unsigned long value)
{
    unsigned char bytes[2];
    bytes[0] = (unsigned char)(value & 0xffUL);
    bytes[1] = (unsigned char)((value >> 8) & 0xffUL);
    fwrite(bytes, 1, 2, file);
}

static void put_u32(FILE *file, unsigned long value)
{
    unsigned char bytes[4];
    bytes[0] = (unsigned char)(value & 0xffUL);
    bytes[1] = (unsigned char)((value >> 8) & 0xffUL);
    bytes[2] = (unsigned char)((value >> 16) & 0xffUL);
    bytes[3] = (unsigned char)((value >> 24) & 0xffUL);
    fwrite(bytes, 1, 4, file);
}

static unsigned long crc32_bytes(const unsigned char *data, size_t size)
{
    unsigned long crc;
    size_t i;
    int bit;

    crc = 0xffffffffUL;
    for (i = 0; i < size; ++i) {
        crc ^= data[i];
        for (bit = 0; bit < 8; ++bit)
            crc = (crc >> 1) ^ (0xedb88320UL & (0UL - (crc & 1UL)));
    }
    return (crc ^ 0xffffffffUL) & 0xffffffffUL;
}

static int host_is_little_endian(void)
{
    const uint16_t one = 1;
    return *((const unsigned char *)&one) == 1;
}

static size_t element_count(int ndim, const size_t *shape)
{
    size_t count;
    int dim;

    count = 1;
    for (dim = 0; dim < ndim; ++dim)
        count *= shape[dim];
    return count;
}

static int make_npy(unsigned char **payload, size_t *payload_size,
                    const char *descr, const void *data, size_t item_size,
                    int ndim, const size_t *shape)
{
    char shape_text[160], dictionary[256], part[40];
    unsigned char *buffer;
    size_t count, data_size, dict_size, header_size, total_size;
    size_t i, j, offset;
    int dim;

    shape_text[0] = '\0';
    strcat(shape_text, "(");
    for (dim = 0; dim < ndim; ++dim) {
        snprintf(part, sizeof(part), "%lu", (unsigned long)shape[dim]);
        strcat(shape_text, part);
        if (dim + 1 < ndim || ndim == 1)
            strcat(shape_text, ",");
        if (dim + 1 < ndim)
            strcat(shape_text, " ");
    }
    strcat(shape_text, ")");

    snprintf(dictionary, sizeof(dictionary),
             "{'descr': '%s', 'fortran_order': False, 'shape': %s, }",
             descr, shape_text);
    dict_size = strlen(dictionary);
    header_size = dict_size + 1;
    while ((10 + header_size) % 16 != 0)
        ++header_size;

    count = element_count(ndim, shape);
    data_size = count * item_size;
    total_size = 10 + header_size + data_size;
    buffer = (unsigned char *)malloc(total_size);
    if (!buffer)
        return -1;

    buffer[0] = 0x93;
    memcpy(buffer + 1, "NUMPY", 5);
    buffer[6] = 1;
    buffer[7] = 0;
    buffer[8] = (unsigned char)(header_size & 0xffUL);
    buffer[9] = (unsigned char)((header_size >> 8) & 0xffUL);
    memcpy(buffer + 10, dictionary, dict_size);
    memset(buffer + 10 + dict_size, ' ', header_size - dict_size);
    buffer[10 + header_size - 1] = '\n';

    offset = 10 + header_size;
    if (host_is_little_endian())
        memcpy(buffer + offset, data, data_size);
    else {
        const unsigned char *source = (const unsigned char *)data;
        for (i = 0; i < count; ++i)
            for (j = 0; j < item_size; ++j)
                buffer[offset + i * item_size + j] =
                    source[i * item_size + item_size - 1 - j];
    }

    *payload = buffer;
    *payload_size = total_size;
    return 0;
}

static int add_array(struct Npz_writer *writer, const char *name,
                     const char *descr, const void *data, size_t item_size,
                     int ndim, const size_t *shape)
{
    struct Npz_entry *entry;
    unsigned char *payload;
    size_t payload_size, filename_size;
    long offset;

    if (!writer->file || writer->failed)
        return -1;
    if (writer->entry_count >= NPZ_MAX_ENTRIES) {
        fprintf(stderr,
                "NPZ entry capacity exceeded while adding '%s': %d >= %d\n",
                name, writer->entry_count, NPZ_MAX_ENTRIES);
        writer->failed = 1;
        return -1;
    }
    if (make_npy(&payload, &payload_size, descr, data, item_size,
                 ndim, shape) != 0) {
        writer->failed = 1;
        return -1;
    }

    entry = &writer->entries[writer->entry_count];
    snprintf(entry->name, sizeof(entry->name), "%s.npy", name);
    entry->crc = crc32_bytes(payload, payload_size);
    entry->size = (unsigned long)payload_size;
    offset = ftell(writer->file);
    if (offset < 0) {
        free(payload);
        writer->failed = 1;
        return -1;
    }
    entry->offset = (unsigned long)offset;
    filename_size = strlen(entry->name);

    put_u32(writer->file, 0x04034b50UL);
    put_u16(writer->file, 20);
    put_u16(writer->file, 0);
    put_u16(writer->file, 0);
    put_u16(writer->file, 0);
    put_u16(writer->file, 0);
    put_u32(writer->file, entry->crc);
    put_u32(writer->file, entry->size);
    put_u32(writer->file, entry->size);
    put_u16(writer->file, filename_size);
    put_u16(writer->file, 0);
    fwrite(entry->name, 1, filename_size, writer->file);
    fwrite(payload, 1, payload_size, writer->file);
    free(payload);

    if (ferror(writer->file)) {
        writer->failed = 1;
        return -1;
    }
    ++writer->entry_count;
    return 0;
}

int npz_open(struct Npz_writer *writer, const char *filename)
{
    memset(writer, 0, sizeof(*writer));
    writer->file = fopen(filename, "wb");
    return writer->file ? 0 : -1;
}

int npz_add_f64(struct Npz_writer *writer, const char *name,
                const double *data, int ndim, const size_t *shape)
{
    return add_array(writer, name, "<f8", data, sizeof(double), ndim, shape);
}

int npz_add_i32(struct Npz_writer *writer, const char *name,
                const int *data, int ndim, const size_t *shape)
{
    return add_array(writer, name, "<i4", data, sizeof(int), ndim, shape);
}

int npz_close(struct Npz_writer *writer)
{
    unsigned long central_offset, central_size;
    long position;
    size_t filename_size;
    int i, failed;

    if (!writer->file)
        return -1;

    position = ftell(writer->file);
    if (position < 0)
        writer->failed = 1;
    central_offset = position < 0 ? 0UL : (unsigned long)position;

    for (i = 0; i < writer->entry_count; ++i) {
        struct Npz_entry *entry = &writer->entries[i];
        filename_size = strlen(entry->name);
        put_u32(writer->file, 0x02014b50UL);
        put_u16(writer->file, 20);
        put_u16(writer->file, 20);
        put_u16(writer->file, 0);
        put_u16(writer->file, 0);
        put_u16(writer->file, 0);
        put_u16(writer->file, 0);
        put_u32(writer->file, entry->crc);
        put_u32(writer->file, entry->size);
        put_u32(writer->file, entry->size);
        put_u16(writer->file, filename_size);
        put_u16(writer->file, 0);
        put_u16(writer->file, 0);
        put_u16(writer->file, 0);
        put_u16(writer->file, 0);
        put_u32(writer->file, 0);
        put_u32(writer->file, entry->offset);
        fwrite(entry->name, 1, filename_size, writer->file);
    }

    position = ftell(writer->file);
    if (position < 0)
        writer->failed = 1;
    central_size = position < 0 ? 0UL :
        (unsigned long)position - central_offset;

    put_u32(writer->file, 0x06054b50UL);
    put_u16(writer->file, 0);
    put_u16(writer->file, 0);
    put_u16(writer->file, writer->entry_count);
    put_u16(writer->file, writer->entry_count);
    put_u32(writer->file, central_size);
    put_u32(writer->file, central_offset);
    put_u16(writer->file, 0);

    failed = writer->failed || ferror(writer->file) || fclose(writer->file) != 0;
    writer->file = NULL;
    return failed ? -1 : 0;
}
