#ifndef CITCOMS_NPZ_WRITER_H
#define CITCOMS_NPZ_WRITER_H

#include <stdio.h>
#include <stddef.h>

#define NPZ_MAX_ENTRIES 512
#define NPZ_MAX_NAME 128

struct Npz_entry {
    char name[NPZ_MAX_NAME];
    unsigned long crc;
    unsigned long size;
    unsigned long offset;
};

struct Npz_writer {
    FILE *file;
    struct Npz_entry entries[NPZ_MAX_ENTRIES];
    int entry_count;
    int failed;
};

int npz_open(struct Npz_writer *, const char *);
int npz_add_f64(struct Npz_writer *, const char *, const double *,
                int, const size_t *);
int npz_add_i32(struct Npz_writer *, const char *, const int *,
                int, const size_t *);
int npz_close(struct Npz_writer *);

#endif
