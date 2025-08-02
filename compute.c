#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FILE_PREFIX "dna-sequences"
#define FILE_EXTENSION ".fasta"
#define ull unsigned long long

struct AminoAcidCodons {
    char amino_acid;
    const char* codons[7];
};

struct AminoAcidCodons codon_map[] = {
    {'F', {"TTT", "TTC", NULL, NULL, NULL, NULL, NULL}},  // Phe
    {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC", NULL}},  // Ser
    {'Y', {"TAT", "TAC", NULL, NULL, NULL, NULL, NULL}},  // Tyr
    {'C', {"TGT", "TGC", NULL, NULL, NULL, NULL, NULL}},  // Cys
    {'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG", NULL}},  // Leu
    {'P', {"CCT", "CCC", "CCA", "CCG", NULL, NULL, NULL}},  // Pro
    {'H', {"CAT", "CAC", NULL, NULL, NULL, NULL, NULL}},  // His
    {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG", NULL}},  // Arg
    {'Q', {"CAA","CAG", NULL, NULL, NULL, NULL, NULL}}, // Gln
    {'I', {"ATT", "ATC", "ATA", NULL, NULL, NULL, NULL}},  // Ile
    {'T', {"ACT", "ACC", "ACA", "ACG", NULL, NULL, NULL}},  // Thr
    {'N', {"AAT", "AAC", NULL, NULL, NULL, NULL, NULL}},  // Asn
    {'K', {"AAA", "AAG", NULL, NULL, NULL, NULL, NULL}},  // Lys
    {'M', {"ATG", NULL, NULL, NULL, NULL, NULL, NULL}},  // Met
    {'V', {"GTT", "GTC", "GTA", "GTG", NULL, NULL, NULL}},  // Val
    {'A', {"GCT", "GCC", "GCA", "GCG", NULL, NULL, NULL}},  // Ala
    {'D', {"GAT", "GAC", NULL, NULL, NULL, NULL, NULL}},  // Asp
    {'G', {"GGT", "GGC", "GGA", "GGG", NULL, NULL, NULL}},  // Gly
    {'E', {"GAA", "GAG", NULL, NULL, NULL, NULL, NULL}},  // Glu
    {'W', {"TGG", NULL, NULL, NULL, NULL, NULL, NULL}},  // Trp
    {'*', {"TAA", "TAG", "TGA", NULL, NULL, NULL, NULL}}  // Stop codons
};

int count_strings(const char** arr) {
    int count = 0;
    while (arr[count] != NULL) {
        count++;
    }
    return count;
}

const char** get_codons(char aminoacid, char** keys, float* values) {
    const char** codons = NULL;

    if (keys == NULL && values == NULL) {
        for (int i = 0; i < sizeof(codon_map) / sizeof(codon_map[0]); i++) {
            if (codon_map[i].amino_acid == aminoacid) {
                return codon_map[i].codons;
            }
        }
        return NULL;  // Return NULL if amino acid is not found
    }

    // Find codons for the amino acid
    for (int i = 0; i < sizeof(codon_map) / sizeof(codon_map[0]); i++) {
        if (codon_map[i].amino_acid == aminoacid) {
            codons = codon_map[i].codons;
            break;
        }
    }

    if (!codons) return NULL;  // Handle missing codons

    // Allocate memory for optimized codon selection
    const char** codons_o = malloc(3 * sizeof(char*));  // Max 2 codons + NULL terminator
    if (!codons_o) return NULL;

    int max1 = 0, max2 = 0;
    int found = 0;

    for (int i = 0; i < count_strings(codons); i++) {
        for (int j = 0; j < count_strings((const char**)keys); j++) {
            if (strcmp(codons[i], keys[j]) == 0) {
                found = 1;
                if (values[j] > values[max1]) {
                    max2 = max1;
                    max1 = j;
                } else if (values[j] > values[max2]) {
                    max2 = j;
                }
            }
        }
    }

    if (!found) {  // If no matching codons were found, return the default codons
        free(codons_o);
        return codons;
    }

    // Assign the best codons
    if (count_strings(codons) == 1) {
        codons_o[0] = keys[max1];
        codons_o[1] = NULL;
    } else {
        codons_o[0] = keys[max1];
        codons_o[1] = keys[max2];
        codons_o[2] = NULL;
    }

    return codons_o;
}


void print_arr(int* arr,int size){
    for (int i=0; i<size; i++){
        printf("%d ",arr[i]);
    }
    printf("\n");
}

__declspec(dllexport) double get_combinations(int optimisation, char* protein) {
    double prod = 0.0;

    for (int i = 0; i < strlen(protein); i++) {
        const char** codons = get_codons(protein[i],NULL,NULL);
        
        int codon_count = count_strings(codons);
        if (optimisation && codon_count > 2) codon_count = 2; 

        prod += log10(codon_count);
    }
    return prod;
}

void generate_sequences(FILE* fptr, char* protein, char* current_seq, int index, int size, int* sequence_count, int max_sequences,char** keys, float* values) {
    if (*sequence_count >= max_sequences) return;

    if (index == size) {
        (*sequence_count)++;
        
        char complementary_seq[size * 3 + 1];
        complementary_seq[0] = '\0';
        for (int i = 0; i < strlen(current_seq); i++) {
            switch (current_seq[i]) {
                case 'A': snprintf(complementary_seq + strlen(complementary_seq), sizeof(complementary_seq) - strlen(complementary_seq), "T"); break;
                case 'T': snprintf(complementary_seq + strlen(complementary_seq), sizeof(complementary_seq) - strlen(complementary_seq), "A"); break;
                case 'C': snprintf(complementary_seq + strlen(complementary_seq), sizeof(complementary_seq) - strlen(complementary_seq), "G"); break;
                case 'G': snprintf(complementary_seq + strlen(complementary_seq), sizeof(complementary_seq) - strlen(complementary_seq), "C"); break;
            }
        }

        char combined_seq[2 * strlen(current_seq) + 1];
        combined_seq[0] = '\0';

        for (int i = 0; i < strlen(current_seq); i++) {
            // Append the base from seq
            snprintf(combined_seq + strlen(combined_seq), sizeof(combined_seq) - strlen(combined_seq), "%c", current_seq[i]);
    
            // Append the corresponding base from complseq
            snprintf(combined_seq + strlen(combined_seq), sizeof(combined_seq) - strlen(combined_seq), "%c", complementary_seq[i]);
        }

        fprintf(fptr, ">Sequence%d\n%s\n", *sequence_count, combined_seq);

        return;
    }

    const char** codons = get_codons(protein[index],keys,values);

    if (!codons) return;

    for (int i = 0; i < count_strings(codons); i++) {
        char new_seq[size * 3 + 1];
        const char* codon = codons[i];
        char* dna = malloc(8);
        dna[0] = '\0';

        for (int x = 0; x < strlen(codon); x++) {
            snprintf(dna, 8, "%s%c", dna, codon[x]);
        }
        
        snprintf(new_seq, sizeof(new_seq), "%s%s", current_seq, dna);  

        generate_sequences(fptr, protein, new_seq, index + 1, size, sequence_count, max_sequences,keys,values);

        free(dna);
    }
}



__declspec(dllexport) void compute(char** keys, float* values, char* protein, int number) {
    int i = 1;
    FILE* file;
    char name[100] = "dna-sequences0.fasta";

    while (1) {  
        if ((file = fopen(name, "r")) != NULL) {  
            fclose(file);
            snprintf(name, sizeof(name), "%s%d%s", FILE_PREFIX, i, FILE_EXTENSION);
            i++;
        } else {
            file = fopen(name, "w");
            if (!file) {
                printf("Error: Cannot create file\n");
                return ;
            }
            break;
        }
    }

    if (keys == NULL && values == NULL) {  
        int sequence_count = 0;
        char current_seq[1] = "";  

        generate_sequences(file, protein, current_seq, 0, strlen(protein), &sequence_count, number,keys,values);
        fclose(file);
    } else {
        int sequence_count = 0;
        char current_seq[1] = "";  

        generate_sequences(file, protein, current_seq, 0, strlen(protein), &sequence_count, number,keys,values);
        fclose(file);
    }
}