#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <libgen.h>
#include <math.h>

#define VERSION "2.alpha"
#define MAX_LINE_LENGTH 65536
#define MAX_PATH_LENGTH 4096
#define MAX_SEQ_LENGTH 100000000
#define DEFAULT_THREADS 8
#define MAX_PRIMER_LENGTH 100
#define BATCH_SIZE 1000

// DNA base bit representation for IUPAC matching
#define A_BIT 0x01
#define C_BIT 0x02
#define G_BIT 0x04
#define T_BIT 0x08

// 2-bit encoding for DNA storage
#define DNA_A 0
#define DNA_C 1
#define DNA_G 2
#define DNA_T 3

// Memory-mapped file structure
typedef struct {
    void *data;
    size_t size;
    int fd;
} MappedFile;

// 2bit file structure
typedef struct {
    uint32_t num_sequences;
    char **headers;
    uint32_t *lengths;
    uint8_t **data;
    uint32_t *packed_lengths;
    MappedFile mapped;
} TwoBitFile;

// Quality file structure (binary format)
typedef struct {
    uint32_t num_sequences;
    uint32_t *lengths;
    uint8_t **data;
    MappedFile mapped;
} QualityFile;

// Primer match structure
typedef struct {
    int type;  // 0=F, 1=F_RC, 2=R, 3=R_RC
    uint32_t start;
    uint32_t end;
    char primer_seq[MAX_PRIMER_LENGTH];
} PrimerMatch;

// Amplicon structure
typedef struct {
    char *sequence;
    char *quality;
    char *source_file;
    char *seq_id;
    uint32_t start;
    uint32_t end;
    char orientation[3];
    char *forward_barcode;
    char *reverse_barcode;
    char forward_primer_seq[MAX_PRIMER_LENGTH];
    char reverse_primer_seq[MAX_PRIMER_LENGTH];
} Amplicon;

// Thread data for compression
typedef struct {
    int thread_id;
    char **files;
    int num_files;
    int start_idx;
    int end_idx;
    char *output_dir;
    FILE *list_fp;
    pthread_mutex_t *list_mutex;
} CompressThreadData;

// Thread data for amplicon finding
typedef struct {
    int thread_id;
    char **seq_files;
    char **qual_files;
    int num_files;
    int start_idx;
    int end_idx;
    
    // Primer info
    uint8_t *forward_mask;
    uint8_t *reverse_mask;
    int forward_len;
    int reverse_len;
    
    // Parameters
    int mismatches;
    int clamp;
    int min_length;
    int max_length;
    int fb_len;
    int rb_len;
    bool trim_primers;
    bool include_offtarget;
    
    // Results
    Amplicon *amplicons;
    int num_amplicons;
    int capacity;
    pthread_mutex_t *result_mutex;
} AmpliconThreadData;

// Global lookup tables
static uint8_t CHAR_TO_MASK[256];
static uint8_t CHAR_TO_2BIT[256];
static char BIT2_TO_BASE[4] = {'A', 'C', 'G', 'T'};

// Initialize lookup tables
void init_lookup_tables() {
    memset(CHAR_TO_MASK, 0, 256);
    memset(CHAR_TO_2BIT, 0, 256);
    
    // Basic nucleotides for bit masks
    CHAR_TO_MASK['A'] = CHAR_TO_MASK['a'] = A_BIT;
    CHAR_TO_MASK['C'] = CHAR_TO_MASK['c'] = C_BIT;
    CHAR_TO_MASK['G'] = CHAR_TO_MASK['g'] = G_BIT;
    CHAR_TO_MASK['T'] = CHAR_TO_MASK['t'] = T_BIT;
    CHAR_TO_MASK['U'] = CHAR_TO_MASK['u'] = T_BIT;
    
    // IUPAC ambiguity codes
    CHAR_TO_MASK['R'] = CHAR_TO_MASK['r'] = A_BIT | G_BIT;
    CHAR_TO_MASK['Y'] = CHAR_TO_MASK['y'] = C_BIT | T_BIT;
    CHAR_TO_MASK['S'] = CHAR_TO_MASK['s'] = G_BIT | C_BIT;
    CHAR_TO_MASK['W'] = CHAR_TO_MASK['w'] = A_BIT | T_BIT;
    CHAR_TO_MASK['K'] = CHAR_TO_MASK['k'] = G_BIT | T_BIT;
    CHAR_TO_MASK['M'] = CHAR_TO_MASK['m'] = A_BIT | C_BIT;
    CHAR_TO_MASK['B'] = CHAR_TO_MASK['b'] = C_BIT | G_BIT | T_BIT;
    CHAR_TO_MASK['D'] = CHAR_TO_MASK['d'] = A_BIT | G_BIT | T_BIT;
    CHAR_TO_MASK['H'] = CHAR_TO_MASK['h'] = A_BIT | C_BIT | T_BIT;
    CHAR_TO_MASK['V'] = CHAR_TO_MASK['v'] = A_BIT | C_BIT | G_BIT;
    CHAR_TO_MASK['N'] = CHAR_TO_MASK['n'] = A_BIT | C_BIT | G_BIT | T_BIT;
    
    // 2-bit encoding
    for (int i = 0; i < 256; i++) {
        CHAR_TO_2BIT[i] = DNA_A;  // Default to A
    }
    CHAR_TO_2BIT['A'] = CHAR_TO_2BIT['a'] = DNA_A;
    CHAR_TO_2BIT['C'] = CHAR_TO_2BIT['c'] = DNA_C;
    CHAR_TO_2BIT['G'] = CHAR_TO_2BIT['g'] = DNA_G;
    CHAR_TO_2BIT['T'] = CHAR_TO_2BIT['t'] = DNA_T;
    CHAR_TO_2BIT['U'] = CHAR_TO_2BIT['u'] = DNA_T;
}

// Convert primer string to bit mask array
uint8_t* primer_to_mask(const char *primer, int *len) {
    *len = strlen(primer);
    uint8_t *mask = (uint8_t*)malloc(*len);
    
    for (int i = 0; i < *len; i++) {
        mask[i] = CHAR_TO_MASK[(unsigned char)primer[i]];
        if (mask[i] == 0) {
            mask[i] = A_BIT | C_BIT | G_BIT | T_BIT;  // Treat unknown as N
        }
    }
    
    return mask;
}

// Reverse complement a bit mask
uint8_t* reverse_complement_mask(const uint8_t *mask, int len) {
    uint8_t *rc = (uint8_t*)malloc(len);
    
    for (int i = 0; i < len; i++) {
        uint8_t orig = mask[len - 1 - i];
        uint8_t comp = 0;
        
        if (orig & A_BIT) comp |= T_BIT;
        if (orig & T_BIT) comp |= A_BIT;
        if (orig & G_BIT) comp |= C_BIT;
        if (orig & C_BIT) comp |= G_BIT;
        
        rc[i] = comp;
    }
    
    return rc;
}

// Reverse complement a DNA string
char* reverse_complement_string(const char *seq, int len) {
    char *rc = (char*)malloc(len + 1);
    
    for (int i = 0; i < len; i++) {
        char c = toupper(seq[len - 1 - i]);
        switch(c) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'G': rc[i] = 'C'; break;
            case 'C': rc[i] = 'G'; break;
            default: rc[i] = 'N'; break;
        }
    }
    rc[len] = '\0';
    
    return rc;
}

// Process spaces in sequence ID
void process_spaces_in_id(char *id) {
    for (int i = 0; id[i]; i++) {
        if (id[i] == ' ' || id[i] == '\t') id[i] = '_';
    }
}

// Check if file has FASTA/FASTQ extension
bool is_sequence_file(const char *filename) {
    const char *ext = strrchr(filename, '.');
    if (!ext) return false;
    
    // Convert to lowercase for comparison
    char lower_ext[10];
    int i;
    for (i = 0; ext[i] && i < 9; i++) {
        lower_ext[i] = tolower(ext[i]);
    }
    lower_ext[i] = '\0';
    
    return (strcmp(lower_ext, ".fa") == 0 ||
            strcmp(lower_ext, ".fasta") == 0 ||
            strcmp(lower_ext, ".fna") == 0 ||
            strcmp(lower_ext, ".ffn") == 0 ||
            strcmp(lower_ext, ".faa") == 0 ||
            strcmp(lower_ext, ".frn") == 0 ||
            strcmp(lower_ext, ".fq") == 0 ||
            strcmp(lower_ext, ".fastq") == 0);
}

// Recursively find all sequence files in directory
void find_sequence_files(const char *dir_path, char ***files, int *num_files, int *capacity) {
    DIR *dir = opendir(dir_path);
    if (!dir) return;
    
    struct dirent *entry;
    struct stat st;
    char full_path[MAX_PATH_LENGTH];
    
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        
        snprintf(full_path, sizeof(full_path), "%s/%s", dir_path, entry->d_name);
        
        if (stat(full_path, &st) == 0) {
            if (S_ISDIR(st.st_mode)) {
                find_sequence_files(full_path, files, num_files, capacity);
            } else if (S_ISREG(st.st_mode) && is_sequence_file(entry->d_name)) {
                if (*num_files >= *capacity) {
                    *capacity *= 2;
                    *files = (char**)realloc(*files, *capacity * sizeof(char*));
                }
                (*files)[(*num_files)++] = strdup(full_path);
            }
        }
    }
    
    closedir(dir);
}

// Compress a FASTA file to 2bit
void compress_fasta(const char *input_file, const char *output_file) {
    FILE *in = fopen(input_file, "r");
    if (!in) return;
    
    FILE *out = fopen(output_file, "wb");
    if (!out) {
        fclose(in);
        return;
    }
    
    // Write placeholder for sequence count
    uint32_t seq_count = 0;
    fwrite(&seq_count, sizeof(uint32_t), 1, out);
    
    char line[MAX_LINE_LENGTH];
    char *header = NULL;
    char *sequence = (char*)malloc(MAX_SEQ_LENGTH);
    uint32_t seq_pos = 0;
    
    while (fgets(line, MAX_LINE_LENGTH, in)) {
        if (line[0] == '>') {
            if (seq_pos > 0 && header) {
                // Process previous sequence
                process_spaces_in_id(header);
                
                // Write header
                uint32_t header_len = strlen(header);
                fwrite(&header_len, sizeof(uint32_t), 1, out);
                fwrite(header, 1, header_len, out);
                
                // Write sequence length
                fwrite(&seq_pos, sizeof(uint32_t), 1, out);
                
                // Pack and write sequence
                uint32_t packed_len = (seq_pos + 3) / 4;
                uint8_t *packed = (uint8_t*)calloc(packed_len, 1);
                
                for (uint32_t i = 0; i < seq_pos; i++) {
                    uint32_t byte_idx = i / 4;
                    uint32_t bit_pos = 6 - 2 * (i % 4);
                    uint8_t value = CHAR_TO_2BIT[(unsigned char)toupper(sequence[i])];
                    packed[byte_idx] |= (value << bit_pos);
                }
                
                fwrite(&packed_len, sizeof(uint32_t), 1, out);
                fwrite(packed, 1, packed_len, out);
                free(packed);
                
                seq_count++;
                seq_pos = 0;
            }
            
            // New header
            if (header) free(header);
            header = strdup(line + 1);
            header[strcspn(header, "\n\r")] = '\0';
        } else {
            // Append to sequence
            int len = strcspn(line, "\n\r");
            if (seq_pos + len < MAX_SEQ_LENGTH) {
                memcpy(sequence + seq_pos, line, len);
                seq_pos += len;
            }
        }
    }
    
    // Process last sequence
    if (seq_pos > 0 && header) {
        process_spaces_in_id(header);
        
        uint32_t header_len = strlen(header);
        fwrite(&header_len, sizeof(uint32_t), 1, out);
        fwrite(header, 1, header_len, out);
        
        fwrite(&seq_pos, sizeof(uint32_t), 1, out);
        
        uint32_t packed_len = (seq_pos + 3) / 4;
        uint8_t *packed = (uint8_t*)calloc(packed_len, 1);
        
        for (uint32_t i = 0; i < seq_pos; i++) {
            uint32_t byte_idx = i / 4;
            uint32_t bit_pos = 6 - 2 * (i % 4);
            uint8_t value = CHAR_TO_2BIT[(unsigned char)toupper(sequence[i])];
            packed[byte_idx] |= (value << bit_pos);
        }
        
        fwrite(&packed_len, sizeof(uint32_t), 1, out);
        fwrite(packed, 1, packed_len, out);
        free(packed);
        
        seq_count++;
    }
    
    // Update sequence count
    fseek(out, 0, SEEK_SET);
    fwrite(&seq_count, sizeof(uint32_t), 1, out);
    
    if (header) free(header);
    free(sequence);
    fclose(in);
    fclose(out);
}

// Compress a FASTQ file to 2bit + quality
void compress_fastq(const char *input_file, const char *seq_output, const char *qual_output) {
    FILE *in = fopen(input_file, "r");
    if (!in) return;
    
    FILE *seq_out = fopen(seq_output, "wb");
    FILE *qual_out = fopen(qual_output, "wb");
    if (!seq_out || !qual_out) {
        fclose(in);
        if (seq_out) fclose(seq_out);
        if (qual_out) fclose(qual_out);
        return;
    }
    
    // Write placeholders
    uint32_t seq_count = 0;
    fwrite(&seq_count, sizeof(uint32_t), 1, seq_out);
    fwrite(&seq_count, sizeof(uint32_t), 1, qual_out);
    
    char line[MAX_LINE_LENGTH];
    char header[MAX_LINE_LENGTH];
    char sequence[MAX_LINE_LENGTH];
    char quality[MAX_LINE_LENGTH];
    int line_num = 0;
    
    while (fgets(line, MAX_LINE_LENGTH, in)) {
        if (line_num % 4 == 0) {  // Header
            strcpy(header, line + 1);
            header[strcspn(header, "\n\r")] = '\0';
            process_spaces_in_id(header);
        } else if (line_num % 4 == 1) {  // Sequence
            strcpy(sequence, line);
            sequence[strcspn(sequence, "\n\r")] = '\0';
        } else if (line_num % 4 == 3) {  // Quality
            strcpy(quality, line);
            quality[strcspn(quality, "\n\r")] = '\0';
            
            uint32_t seq_len = strlen(sequence);
            
            // Write to 2bit file
            uint32_t header_len = strlen(header);
            fwrite(&header_len, sizeof(uint32_t), 1, seq_out);
            fwrite(header, 1, header_len, seq_out);
            
            fwrite(&seq_len, sizeof(uint32_t), 1, seq_out);
            
            uint32_t packed_len = (seq_len + 3) / 4;
            uint8_t *packed = (uint8_t*)calloc(packed_len, 1);
            
            for (uint32_t i = 0; i < seq_len; i++) {
                uint32_t byte_idx = i / 4;
                uint32_t bit_pos = 6 - 2 * (i % 4);
                uint8_t value = CHAR_TO_2BIT[(unsigned char)toupper(sequence[i])];
                packed[byte_idx] |= (value << bit_pos);
            }
            
            fwrite(&packed_len, sizeof(uint32_t), 1, seq_out);
            fwrite(packed, 1, packed_len, seq_out);
            free(packed);
            
            // Write to quality file
            fwrite(&seq_len, sizeof(uint32_t), 1, qual_out);
            fwrite(quality, 1, seq_len, qual_out);
            
            seq_count++;
        }
        line_num++;
    }
    
    // Update counts
    fseek(seq_out, 0, SEEK_SET);
    fwrite(&seq_count, sizeof(uint32_t), 1, seq_out);
    fseek(qual_out, 0, SEEK_SET);
    fwrite(&seq_count, sizeof(uint32_t), 1, qual_out);
    
    fclose(in);
    fclose(seq_out);
    fclose(qual_out);
}

// Worker thread for compression
void* compress_worker(void *arg) {
    CompressThreadData *data = (CompressThreadData*)arg;
    
    for (int i = data->start_idx; i < data->end_idx; i++) {
        char *input_file = data->files[i];
        
        // Determine if FASTA or FASTQ
        FILE *test = fopen(input_file, "r");
        if (!test) continue;
        
        char first_char = fgetc(test);
        fclose(test);
        
        // Generate output filename
        char *base_name = basename(input_file);
        char seq_output[MAX_PATH_LENGTH];
        char qual_output[MAX_PATH_LENGTH];
        
        // Remove extension
        char name_no_ext[MAX_PATH_LENGTH];
        strcpy(name_no_ext, base_name);
        char *dot = strrchr(name_no_ext, '.');
        if (dot) *dot = '\0';
        
        snprintf(seq_output, sizeof(seq_output), "%s/%s.2bit", data->output_dir, name_no_ext);
        snprintf(qual_output, sizeof(qual_output), "%s/%s.qual", data->output_dir, name_no_ext);
        
        if (first_char == '>') {
            // FASTA
            compress_fasta(input_file, seq_output);
            
            // Write to file list (single column)
            pthread_mutex_lock(data->list_mutex);
            fprintf(data->list_fp, "%s\n", seq_output);
            pthread_mutex_unlock(data->list_mutex);
            
        } else if (first_char == '@') {
            // FASTQ
            compress_fastq(input_file, seq_output, qual_output);
            
            // Write to file list (two columns)
            pthread_mutex_lock(data->list_mutex);
            fprintf(data->list_fp, "%s\t%s\n", seq_output, qual_output);
            pthread_mutex_unlock(data->list_mutex);
        }
    }
    
    return NULL;
}

// Memory map a file
MappedFile map_file(const char *filename) {
    MappedFile mapped = {0};
    mapped.fd = open(filename, O_RDONLY);
    if (mapped.fd == -1) return mapped;
    
    struct stat sb;
    if (fstat(mapped.fd, &sb) == -1) {
        close(mapped.fd);
        mapped.fd = -1;
        return mapped;
    }
    
    mapped.size = sb.st_size;
    mapped.data = mmap(NULL, mapped.size, PROT_READ, MAP_PRIVATE, mapped.fd, 0);
    
    if (mapped.data == MAP_FAILED) {
        close(mapped.fd);
        mapped.fd = -1;
        mapped.data = NULL;
    }
    
    return mapped;
}

// Unmap file
void unmap_file(MappedFile *mapped) {
    if (mapped->data && mapped->data != MAP_FAILED) {
        munmap(mapped->data, mapped->size);
    }
    if (mapped->fd != -1) {
        close(mapped->fd);
    }
}

// Load 2bit file
TwoBitFile load_2bit_file(const char *filename) {
    TwoBitFile file = {0};
    file.mapped = map_file(filename);
    if (!file.mapped.data) return file;
    
    const uint8_t *data = (const uint8_t*)file.mapped.data;
    size_t offset = 0;
    
    file.num_sequences = *(uint32_t*)(data + offset);
    offset += sizeof(uint32_t);
    
    file.headers = (char**)malloc(file.num_sequences * sizeof(char*));
    file.lengths = (uint32_t*)malloc(file.num_sequences * sizeof(uint32_t));
    file.data = (uint8_t**)malloc(file.num_sequences * sizeof(uint8_t*));
    file.packed_lengths = (uint32_t*)malloc(file.num_sequences * sizeof(uint32_t));
    
    for (uint32_t i = 0; i < file.num_sequences; i++) {
        uint32_t header_len = *(uint32_t*)(data + offset);
        offset += sizeof(uint32_t);
        
        file.headers[i] = (char*)malloc(header_len + 1);
        memcpy(file.headers[i], data + offset, header_len);
        file.headers[i][header_len] = '\0';
        offset += header_len;
        
        file.lengths[i] = *(uint32_t*)(data + offset);
        offset += sizeof(uint32_t);
        
        file.packed_lengths[i] = *(uint32_t*)(data + offset);
        offset += sizeof(uint32_t);
        
        file.data[i] = (uint8_t*)(data + offset);
        offset += file.packed_lengths[i];
    }
    
    return file;
}

// Load quality file
QualityFile load_quality_file(const char *filename) {
    QualityFile file = {0};
    file.mapped = map_file(filename);
    if (!file.mapped.data) return file;
    
    const uint8_t *data = (const uint8_t*)file.mapped.data;
    size_t offset = 0;
    
    file.num_sequences = *(uint32_t*)(data + offset);
    offset += sizeof(uint32_t);
    
    file.lengths = (uint32_t*)malloc(file.num_sequences * sizeof(uint32_t));
    file.data = (uint8_t**)malloc(file.num_sequences * sizeof(uint8_t*));
    
    for (uint32_t i = 0; i < file.num_sequences; i++) {
        file.lengths[i] = *(uint32_t*)(data + offset);
        offset += sizeof(uint32_t);
        
        file.data[i] = (uint8_t*)(data + offset);
        offset += file.lengths[i];
    }
    
    return file;
}

// Free 2bit file
void free_2bit_file(TwoBitFile *file) {
    for (uint32_t i = 0; i < file->num_sequences; i++) {
        free(file->headers[i]);
    }
    free(file->headers);
    free(file->lengths);
    free(file->data);
    free(file->packed_lengths);
    unmap_file(&file->mapped);
}

// Free quality file
void free_quality_file(QualityFile *file) {
    free(file->lengths);
    free(file->data);
    unmap_file(&file->mapped);
}

// Unpack 2bit sequence
void unpack_sequence(const uint8_t *packed, uint32_t seq_length, char *output) {
    for (uint32_t i = 0; i < seq_length; i++) {
        uint32_t byte_idx = i / 4;
        uint32_t bit_pos = 6 - 2 * (i % 4);
        uint8_t value = (packed[byte_idx] >> bit_pos) & 0x3;
        output[i] = BIT2_TO_BASE[value];
    }
    output[seq_length] = '\0';
}

// Check if sequence matches primer with mismatches
bool matches_with_mismatches(const char *seq, const uint8_t *primer_mask, 
                             int len, int max_mismatches, int clamp, bool check_rc_clamp) {
    int mismatches = 0;
    
    for (int i = 0; i < len; i++) {
        uint8_t seq_mask = CHAR_TO_MASK[(unsigned char)seq[i]];
        
        if ((primer_mask[i] & seq_mask) == 0) {
            mismatches++;
            
            if (check_rc_clamp) {
                if (i < clamp) return false;
            } else {
                if (i >= len - clamp) return false;
            }
            
            if (mismatches > max_mismatches) return false;
        }
    }
    
    return true;
}

// Find all primer matches in a sequence
PrimerMatch* find_primer_matches(const char *sequence, int seq_len,
                                 const uint8_t *forward_mask, int f_len,
                                 const uint8_t *reverse_mask, int r_len,
                                 int max_mismatches, int clamp, int *num_matches) {
    uint8_t *forward_rc_mask = reverse_complement_mask(forward_mask, f_len);
    uint8_t *reverse_rc_mask = reverse_complement_mask(reverse_mask, r_len);
    
    int capacity = 100;
    PrimerMatch *matches = (PrimerMatch*)malloc(capacity * sizeof(PrimerMatch));
    *num_matches = 0;
    
    // Search for forward primer
    for (int pos = 0; pos <= seq_len - f_len; pos++) {
        if (matches_with_mismatches(sequence + pos, forward_mask, f_len, 
                                   max_mismatches, clamp, false)) {
            if (*num_matches >= capacity) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 0;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + f_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, f_len);
            matches[*num_matches].primer_seq[f_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Search for forward primer RC
    for (int pos = 0; pos <= seq_len - f_len; pos++) {
        if (matches_with_mismatches(sequence + pos, forward_rc_mask, f_len, 
                                   max_mismatches, clamp, true)) {
            if (*num_matches >= capacity) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 1;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + f_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, f_len);
            matches[*num_matches].primer_seq[f_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Search for reverse primer
    for (int pos = 0; pos <= seq_len - r_len; pos++) {
        if (matches_with_mismatches(sequence + pos, reverse_mask, r_len, 
                                   max_mismatches, clamp, false)) {
            if (*num_matches >= capacity) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 2;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + r_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, r_len);
            matches[*num_matches].primer_seq[r_len] = '\0';
            (*num_matches)++;
        }
    }
    
    // Search for reverse primer RC
    for (int pos = 0; pos <= seq_len - r_len; pos++) {
        if (matches_with_mismatches(sequence + pos, reverse_rc_mask, r_len, 
                                   max_mismatches, clamp, true)) {
            if (*num_matches >= capacity) {
                capacity *= 2;
                matches = (PrimerMatch*)realloc(matches, capacity * sizeof(PrimerMatch));
            }
            
            matches[*num_matches].type = 3;
            matches[*num_matches].start = pos;
            matches[*num_matches].end = pos + r_len;
            strncpy(matches[*num_matches].primer_seq, sequence + pos, r_len);
            matches[*num_matches].primer_seq[r_len] = '\0';
            (*num_matches)++;
        }
    }
    
    free(forward_rc_mask);
    free(reverse_rc_mask);
    
    return matches;
}

// Process a sequence to find amplicons
void process_sequence_for_amplicons(const char *seq_id, const char *sequence, const char *quality,
                                    const char *source_file, AmpliconThreadData *data) {
    int seq_len = strlen(sequence);
    
    int num_matches;
    PrimerMatch *matches = find_primer_matches(sequence, seq_len,
                                               data->forward_mask, data->forward_len,
                                               data->reverse_mask, data->reverse_len,
                                               data->mismatches, data->clamp, &num_matches);
    
    // Sort matches by position
    for (int i = 0; i < num_matches - 1; i++) {
        for (int j = i + 1; j < num_matches; j++) {
            if (matches[j].start < matches[i].start) {
                PrimerMatch temp = matches[i];
                matches[i] = matches[j];
                matches[j] = temp;
            }
        }
    }
    
    // Find valid amplicons
    for (int i = 0; i < num_matches - 1; i++) {
        if (matches[i].type != 0 && matches[i].type != 2) continue;
        
        for (int j = i + 1; j < num_matches; j++) {
            int length = matches[j].end - matches[i].start;
            
            if (length > data->max_length) break;
            if (length < data->min_length) continue;
            
            if (matches[j].type != 1 && matches[j].type != 3) continue;
            
            char orientation[3];
            if (matches[i].type == 0 && matches[j].type == 3) {
                strcpy(orientation, "FR");
            } else if (matches[i].type == 2 && matches[j].type == 1) {
                strcpy(orientation, "RF");
            } else if (matches[i].type == 0 && matches[j].type == 1) {
                strcpy(orientation, "FF");
            } else if (matches[i].type == 2 && matches[j].type == 3) {
                strcpy(orientation, "RR");
            } else {
                continue;
            }
            
            if (!data->include_offtarget && strcmp(orientation, "FR") != 0 && strcmp(orientation, "RF") != 0) {
                continue;
            }
            
            pthread_mutex_lock(data->result_mutex);
            
            if (data->num_amplicons >= data->capacity) {
                data->capacity *= 2;
                data->amplicons = (Amplicon*)realloc(data->amplicons, data->capacity * sizeof(Amplicon));
            }
            
            Amplicon *amp = &data->amplicons[data->num_amplicons];
            
            int amp_start = data->trim_primers ? matches[i].end : matches[i].start;
            int amp_end = data->trim_primers ? matches[j].start : matches[j].end;
            int amp_len = amp_end - amp_start;
            
            if (amp_len <= 0) {
                pthread_mutex_unlock(data->result_mutex);
                continue;
            }
            
            amp->sequence = (char*)malloc(amp_len + 1);
            strncpy(amp->sequence, sequence + amp_start, amp_len);
            amp->sequence[amp_len] = '\0';
            
            if (quality) {
                amp->quality = (char*)malloc(amp_len + 1);
                strncpy(amp->quality, quality + amp_start, amp_len);
                amp->quality[amp_len] = '\0';
            } else {
                amp->quality = NULL;
            }
            
            if (strcmp(orientation, "RF") == 0) {
                char *rc = reverse_complement_string(amp->sequence, amp_len);
                free(amp->sequence);
                amp->sequence = rc;
                
                if (amp->quality) {
                    for (int k = 0; k < amp_len / 2; k++) {
                        char temp = amp->quality[k];
                        amp->quality[k] = amp->quality[amp_len - 1 - k];
                        amp->quality[amp_len - 1 - k] = temp;
                    }
                }
            }
            
            // Extract barcodes
            amp->forward_barcode = NULL;
            amp->reverse_barcode = NULL;
            
            if (strcmp(orientation, "FR") == 0) {
                if (data->fb_len > 0 && matches[i].start >= data->fb_len) {
                    amp->forward_barcode = (char*)malloc(data->fb_len + 1);
                    strncpy(amp->forward_barcode, sequence + matches[i].start - data->fb_len, data->fb_len);
                    amp->forward_barcode[data->fb_len] = '\0';
                }
                if (data->rb_len > 0 && matches[j].end + data->rb_len <= seq_len) {
                    amp->reverse_barcode = (char*)malloc(data->rb_len + 1);
                    strncpy(amp->reverse_barcode, sequence + matches[j].end, data->rb_len);
                    amp->reverse_barcode[data->rb_len] = '\0';
                }
            } else if (strcmp(orientation, "RF") == 0) {
                if (data->fb_len > 0 && matches[j].end + data->fb_len <= seq_len) {
                    char *temp = (char*)malloc(data->fb_len + 1);
                    strncpy(temp, sequence + matches[j].end, data->fb_len);
                    temp[data->fb_len] = '\0';
                    amp->forward_barcode = reverse_complement_string(temp, data->fb_len);
                    free(temp);
                }
                if (data->rb_len > 0 && matches[i].start >= data->rb_len) {
                    char *temp = (char*)malloc(data->rb_len + 1);
                    strncpy(temp, sequence + matches[i].start - data->rb_len, data->rb_len);
                    temp[data->rb_len] = '\0';
                    amp->reverse_barcode = reverse_complement_string(temp, data->rb_len);
                    free(temp);
                }
            }
            
            amp->source_file = strdup(source_file);
            amp->seq_id = strdup(seq_id);
            amp->start = matches[i].start;
            amp->end = matches[j].end;
            strcpy(amp->orientation, orientation);
            strcpy(amp->forward_primer_seq, matches[i].primer_seq);
            strcpy(amp->reverse_primer_seq, matches[j].primer_seq);
            
            data->num_amplicons++;
            
            pthread_mutex_unlock(data->result_mutex);
        }
    }
    
    free(matches);
}

// Worker thread for amplicon finding
void* amplicon_worker(void *arg) {
    AmpliconThreadData *data = (AmpliconThreadData*)arg;
    
    for (int idx = data->start_idx; idx < data->end_idx; idx++) {
        TwoBitFile seq_file = load_2bit_file(data->seq_files[idx]);
        if (!seq_file.mapped.data) continue;
        
        QualityFile qual_file = {0};
        if (data->qual_files && data->qual_files[idx]) {
            qual_file = load_quality_file(data->qual_files[idx]);
        }
        
        char *base_name = basename(data->seq_files[idx]);
        
        for (uint32_t i = 0; i < seq_file.num_sequences; i++) {
            char *sequence = (char*)malloc(seq_file.lengths[i] + 1);
            unpack_sequence(seq_file.data[i], seq_file.lengths[i], sequence);
            
            char *quality = NULL;
            if (qual_file.data && i < qual_file.num_sequences) {
                quality = (char*)malloc(qual_file.lengths[i] + 1);
                memcpy(quality, qual_file.data[i], qual_file.lengths[i]);
                quality[qual_file.lengths[i]] = '\0';
            }
            
            process_sequence_for_amplicons(seq_file.headers[i], sequence, quality, base_name, data);
            
            free(sequence);
            if (quality) free(quality);
        }
        
        free_2bit_file(&seq_file);
        if (qual_file.data) free_quality_file(&qual_file);
    }
    
    return NULL;
}

// Compress command
int compress_command(const char *input_dir, const char *output_dir, int num_threads) {
    // Create output directory
    mkdir(output_dir, 0755);
    
    // Find all sequence files
    int capacity = 1000;
    char **files = (char**)malloc(capacity * sizeof(char*));
    int num_files = 0;
    
    fprintf(stderr, "Scanning directory %s for sequence files...\n", input_dir);
    find_sequence_files(input_dir, &files, &num_files, &capacity);
    fprintf(stderr, "Found %d sequence files\n", num_files);
    
    if (num_files == 0) {
        fprintf(stderr, "No sequence files found\n");
        free(files);
        return 1;
    }
    
    // Create temporary file for list
    char list_file[MAX_PATH_LENGTH];
    snprintf(list_file, sizeof(list_file), "%s/file_list.tmp", output_dir);
    FILE *list_fp = fopen(list_file, "w");
    if (!list_fp) {
        fprintf(stderr, "Cannot create file list\n");
        return 1;
    }
    
    pthread_mutex_t list_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    // Create threads
    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    CompressThreadData *thread_data = (CompressThreadData*)malloc(num_threads * sizeof(CompressThreadData));
    
    int files_per_thread = (num_files + num_threads - 1) / num_threads;
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].files = files;
        thread_data[i].num_files = num_files;
        thread_data[i].start_idx = i * files_per_thread;
        thread_data[i].end_idx = ((i + 1) * files_per_thread < num_files) ? 
                                  (i + 1) * files_per_thread : num_files;
        thread_data[i].output_dir = (char*)output_dir;
        thread_data[i].list_fp = list_fp;
        thread_data[i].list_mutex = &list_mutex;
        
        pthread_create(&threads[i], NULL, compress_worker, &thread_data[i]);
    }
    
    // Wait for threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    fclose(list_fp);
    
    // Output file list to stdout
    FILE *read_list = fopen(list_file, "r");
    if (read_list) {
        char line[MAX_PATH_LENGTH * 2];
        while (fgets(line, sizeof(line), read_list)) {
            printf("%s", line);
        }
        fclose(read_list);
        remove(list_file);
    }
    
    // Clean up
    for (int i = 0; i < num_files; i++) {
        free(files[i]);
    }
    free(files);
    free(threads);
    free(thread_data);
    pthread_mutex_destroy(&list_mutex);
    
    fprintf(stderr, "Compression complete. Compressed %d files\n", num_files);
    
    return 0;
}

// Run command
// Replace the run_command function starting at line 1033 with this fixed version:

int run_command(const char *input_list, const char *primers_file, const char *output_file,
                int num_threads, int mismatches, int clamp, int min_length, int max_length,
                int fb_len, int rb_len, bool trim_primers, bool include_offtarget) {
    
    FILE *fp = fopen(primers_file, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open primers file %s\n", primers_file);
        return 1;
    }
    
    char forward_primer[256], reverse_primer[256];
    if (fscanf(fp, "%s %s", forward_primer, reverse_primer) != 2) {
        fprintf(stderr, "Error reading primers from file\n");
        fclose(fp);
        return 1;
    }
    fclose(fp);
    
    fprintf(stderr, "Forward primer: %s\n", forward_primer);
    fprintf(stderr, "Reverse primer: %s\n", reverse_primer);
    
    int forward_len, reverse_len;
    uint8_t *forward_mask = primer_to_mask(forward_primer, &forward_len);
    uint8_t *reverse_mask = primer_to_mask(reverse_primer, &reverse_len);
    
    int capacity = 1000;
    char **seq_files = (char**)malloc(capacity * sizeof(char*));
    char **qual_files = (char**)malloc(capacity * sizeof(char*));
    int num_files = 0;
    
    fp = fopen(input_list, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open input list %s\n", input_list);
        return 1;
    }
    
    char line[MAX_PATH_LENGTH * 2];
    while (fgets(line, sizeof(line), fp)) {
        if (num_files >= capacity) {
            capacity *= 2;
            seq_files = (char**)realloc(seq_files, capacity * sizeof(char*));
            qual_files = (char**)realloc(qual_files, capacity * sizeof(char*));
        }
        
        char *tab = strchr(line, '\t');
        if (tab) {
            *tab = '\0';
            seq_files[num_files] = strdup(line);
            tab++;
            tab[strcspn(tab, "\n\r")] = '\0';
            qual_files[num_files] = strdup(tab);
        } else {
            line[strcspn(line, "\n\r")] = '\0';
            seq_files[num_files] = strdup(line);
            qual_files[num_files] = NULL;
        }
        num_files++;
    }
    fclose(fp);
    
    fprintf(stderr, "Processing %d files with %d threads\n", num_files, num_threads);
    
    pthread_mutex_t result_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    AmpliconThreadData *thread_data = (AmpliconThreadData*)malloc(num_threads * sizeof(AmpliconThreadData));
    
    int files_per_thread = (num_files + num_threads - 1) / num_threads;
    
    // CRITICAL FIX: Give each thread its own amplicons array
    // and merge them at the end
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].seq_files = seq_files;
        thread_data[i].qual_files = qual_files;
        thread_data[i].num_files = num_files;
        thread_data[i].start_idx = i * files_per_thread;
        thread_data[i].end_idx = ((i + 1) * files_per_thread < num_files) ? 
                                  (i + 1) * files_per_thread : num_files;
        
        thread_data[i].forward_mask = forward_mask;
        thread_data[i].reverse_mask = reverse_mask;
        thread_data[i].forward_len = forward_len;
        thread_data[i].reverse_len = reverse_len;
        
        thread_data[i].mismatches = mismatches;
        thread_data[i].clamp = clamp;
        thread_data[i].min_length = min_length;
        thread_data[i].max_length = max_length;
        thread_data[i].fb_len = fb_len;
        thread_data[i].rb_len = rb_len;
        thread_data[i].trim_primers = trim_primers;
        thread_data[i].include_offtarget = include_offtarget;
        
        // CRITICAL FIX: Each thread gets its own amplicons array
        thread_data[i].amplicons = (Amplicon*)malloc(10000 * sizeof(Amplicon));
        thread_data[i].num_amplicons = 0;
        thread_data[i].capacity = 10000;
        thread_data[i].result_mutex = &result_mutex;
        
        pthread_create(&threads[i], NULL, amplicon_worker, &thread_data[i]);
    }
    
    // Wait for threads
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    // Merge all thread results
    int total_amplicons = 0;
    for (int i = 0; i < num_threads; i++) {
        total_amplicons += thread_data[i].num_amplicons;
    }
    
    fprintf(stderr, "Found %d amplicons\n", total_amplicons);
    
    FILE *out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "Cannot create output file %s\n", output_file);
        return 1;
    }
    
    bool is_fastq = false;
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < thread_data[t].num_amplicons; i++) {
            if (thread_data[t].amplicons[i].quality) {
                is_fastq = true;
                break;
            }
        }
        if (is_fastq) break;
    }
    
    // Write amplicons from all threads
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < thread_data[t].num_amplicons; i++) {
            Amplicon *amp = &thread_data[t].amplicons[i];
            
            if (is_fastq) {
                fprintf(out, "@%s.source=%s.coordinates=%u-%u.orientation=%s",
                       amp->seq_id, amp->source_file, amp->start, amp->end, amp->orientation);
            } else {
                fprintf(out, ">%s.source=%s.coordinates=%u-%u.orientation=%s",
                       amp->seq_id, amp->source_file, amp->start, amp->end, amp->orientation);
            }
            
            fprintf(out, ".fprimer=%s.rprimer=%s", amp->forward_primer_seq, amp->reverse_primer_seq);
            
            if (amp->forward_barcode) {
                fprintf(out, ".fb=%s", amp->forward_barcode);
            }
            if (amp->reverse_barcode) {
                fprintf(out, ".rb=%s", amp->reverse_barcode);
            }
            fprintf(out, "\n");
            
            fprintf(out, "%s\n", amp->sequence);
            
            if (is_fastq && amp->quality) {
                fprintf(out, "+\n%s\n", amp->quality);
            }
        }
    }
    
    fclose(out);
    
    // Clean up - each thread has its own amplicons array now
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < thread_data[t].num_amplicons; i++) {
            free(thread_data[t].amplicons[i].sequence);
            if (thread_data[t].amplicons[i].quality) free(thread_data[t].amplicons[i].quality);
            free(thread_data[t].amplicons[i].source_file);
            free(thread_data[t].amplicons[i].seq_id);
            if (thread_data[t].amplicons[i].forward_barcode) free(thread_data[t].amplicons[i].forward_barcode);
            if (thread_data[t].amplicons[i].reverse_barcode) free(thread_data[t].amplicons[i].reverse_barcode);
        }
        free(thread_data[t].amplicons);
    }
    
    free(forward_mask);
    free(reverse_mask);
    
    for (int i = 0; i < num_files; i++) {
        free(seq_files[i]);
        if (qual_files[i]) free(qual_files[i]);
    }
    free(seq_files);
    free(qual_files);
    free(threads);
    free(thread_data);
    pthread_mutex_destroy(&result_mutex);
    
    return 0;
}

void print_usage() {
    printf("AmpliconHunter v%s\n\n", VERSION);
    printf("Usage:\n");
    printf("  amplicon_hunter compress --input-dir <dir> --output <dir> [--threads <n>]\n");
    printf("  amplicon_hunter run --input <file_list> --primers <file> --output <file> [options]\n\n");
    
    printf("Compress options:\n");
    printf("  --input-dir <dir>    Directory containing FASTA/FASTQ files\n");
    printf("  --output <dir>       Output directory for compressed files\n");
    printf("  --threads <n>        Number of threads (default: 8)\n\n");
    
    printf("Run options:\n");
    printf("  --input <file>       File list from compress (1 or 2 columns)\n");
    printf("  --primers <file>     Primers file (forward reverse)\n");
    printf("  --output <file>      Output FASTA/FASTQ file\n");
    printf("  --threads <n>        Number of threads (default: 8)\n");
    printf("  --mismatches <n>     Maximum mismatches (default: 0)\n");
    printf("  --clamp <n>          3' clamp size (default: 5)\n");
    printf("  --min-length <n>     Minimum amplicon length (default: 50)\n");
    printf("  --max-length <n>     Maximum amplicon length (default: 5000)\n");
    printf("  --fb-len <n>         Forward barcode length (default: 0)\n");
    printf("  --rb-len <n>         Reverse barcode length (default: 0)\n");
    printf("  --trim-primers       Trim primer sequences\n");
    printf("  --include-offtarget  Include FF and RR orientations\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    init_lookup_tables();
    
    if (strcmp(argv[1], "compress") == 0) {
        char *input_dir = NULL;
        char *output_dir = NULL;
        int num_threads = DEFAULT_THREADS;
        
        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "--input-dir") == 0 && i + 1 < argc) {
                input_dir = argv[++i];
            } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
                output_dir = argv[++i];
            } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
                num_threads = atoi(argv[++i]);
            }
        }
        
        if (!input_dir || !output_dir) {
            fprintf(stderr, "Error: Missing required arguments\n");
            print_usage();
            return 1;
        }
        
        return compress_command(input_dir, output_dir, num_threads);
        
    } else if (strcmp(argv[1], "run") == 0) {
        char *input_list = NULL;
        char *primers_file = NULL;
        char *output_file = NULL;
        int num_threads = DEFAULT_THREADS;
        int mismatches = 0;
        int clamp = 5;
        int min_length = 50;
        int max_length = 5000;
        int fb_len = 0;
        int rb_len = 0;
        bool trim_primers = false;
        bool include_offtarget = false;
        
        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "--input") == 0 && i + 1 < argc) {
                input_list = argv[++i];
            } else if (strcmp(argv[i], "--primers") == 0 && i + 1 < argc) {
                primers_file = argv[++i];
            } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
                output_file = argv[++i];
            } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
                num_threads = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--mismatches") == 0 && i + 1 < argc) {
                mismatches = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--clamp") == 0 && i + 1 < argc) {
                clamp = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--min-length") == 0 && i + 1 < argc) {
                min_length = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--max-length") == 0 && i + 1 < argc) {
                max_length = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--fb-len") == 0 && i + 1 < argc) {
                fb_len = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--rb-len") == 0 && i + 1 < argc) {
                rb_len = atoi(argv[++i]);
            } else if (strcmp(argv[i], "--trim-primers") == 0) {
                trim_primers = true;
            } else if (strcmp(argv[i], "--include-offtarget") == 0) {
                include_offtarget = true;
            }
        }
        
        if (!input_list || !primers_file || !output_file) {
            fprintf(stderr, "Error: Missing required arguments\n");
            print_usage();
            return 1;
        }
        
        return run_command(input_list, primers_file, output_file, num_threads,
                          mismatches, clamp, min_length, max_length,
                          fb_len, rb_len, trim_primers, include_offtarget);
        
    } else {
        fprintf(stderr, "Error: Unknown command '%s'\n", argv[1]);
        print_usage();
        return 1;
    }
    
    return 0;
}