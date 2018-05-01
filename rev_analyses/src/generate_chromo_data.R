
#
# Generate chromosome data input file for RevBayes analysis.
#
# Will Freyman
#

in_file = "../csomeseries/soldataalltaxa.csv"
out_file = "data/chromo_counts.tsv"
out_file_notes = "data/chromo_counts_notes.txt"

d = read.csv(in_file, stringsAsFactors=FALSE)
a = FALSE
a_notes = FALSE
max_c = 0
min_c = 1000

for (i in 1:nrow(d)) {

    counts = vector()

    # chromo counts are in columns 2 thru 16
    # of original data file
    for (j in 2:16) {
        
        c = d[i,j]
        if (!is.na(c) & c %% 2 == 0) {
            
            c = c/2
            counts = c(counts, c)
            if (c > max_c)
                max_c = c
            if (c < min_c)
                min_c = c

        } else if (!is.na(c) & c %% 2 != 0) {
    
            write(paste0(d[i,1], ": dropped odd diploid chromo count = ", c),
                  file=out_file_notes,
                  append=a_notes)
            a_notes = TRUE

        }

    }

    if (length(counts) == 0) {

        counts = paste0(d[i,1], "\t", "?")

    } else if (length(counts) == 1) {
        
        counts = paste0(d[i,1], "\t", counts)

    } else {

        counts = paste0(counts, collapse=" ")
        counts = paste0(d[i,1], "\t", "(", counts, ")")

    }

    write(counts, file=out_file, append=a)
    a = TRUE

}

write(paste0("min haploid chromosome count = ", min_c),
             file=out_file_notes,
             append=a_notes)
write(paste0("max haploid chromosome count = ", max_c),
             file=out_file_notes,
             append=a_notes)
