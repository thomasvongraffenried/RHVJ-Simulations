# Load libraries

library(ggplot2)
library(gridExtra)
library(plotrix)
library(ggforce)
library(ape)
library(geosphere)
library(vegan)
library(dplyr)
library(patchwork)
library(boot)
library(seqinr)
library(adegenet)
library(stringr)
library(caroline)
library(igraph)
library(future)
library(future.apply)
library(pbapply)

# Parameters
n_generations <- 300
space_size <- 1
circle_radius <- 0.001
mutation_rate <- 0.001
K <- 1000
init_seq <- "CCGTATATGCAGAAGGCATACAACATCGCGCCTAGCGTGTATACGGGCGAGACATCTATTGGCACTGGCACGAAACTGACGTACGCAACGTACGGCAAGGCTGCCGCAATGGACACTAGCCTCCTCACGGGCTATGATGTTGTCATCTGTGATGAGTGCCATGACGTTACTGCCACCACCATTCTGGGAATTGGACACGTGCTGACGAAGGCTGAGTCCTGTGGGGTTAAGCTTGTCATTCTGGCCACGGCAACGCCTCCCGGATGTTCAACAACACCACATCCGAACATCACGGAGGTCGAGCTAGGCTCGAGTGGCGAGGTGCAATTCTACGGGAAGCGACTGGAGTTGGCACACTACCTAAAAGGCAGGCATCTCATCTTCTGCGCGTCAAAGCTCGTCTGTGACACCTTGGCGAGCTTGTTGCGGCAACATGGCATTACTGCAGTGGCCTACTACCGAGGTGAA"
nucleotides <- c("A", "C", "G", "T")
seeds <- 1:27  

# Function to mutate sequences
mutate_sequences <- function(seqs, mu) {
  seq_len <- nchar(seqs[1])
  seq_matrix <- do.call(rbind, strsplit(seqs, ""))
  mutation_mask <- matrix(runif(nrow(seq_matrix) * seq_len) < mu, nrow = nrow(seq_matrix))
  seq_matrix[mutation_mask] <- sample(nucleotides, sum(mutation_mask), replace = TRUE)
  apply(seq_matrix, 1, paste, collapse = "")
}

# Loop over seeds
for (seed in seeds) {
  set.seed(seed)
  start_year <- 2025 - n_generations
  cat("Running simulation for seed:", seed, "\n")
  
  # Initial individual
  individuals <- data.frame(
    ID = 1, ParentID = NA, Generation = 1,
    x = 0.5, y = 0.5, Sequence = init_seq, stringsAsFactors = FALSE
  )
  individuals$Sequence_name <- paste0(individuals$ID, "_", start_year)
  all_individuals_list <- list(individuals)
  next_ID <- 2
  
  # Simulation loop
  for (gen in 2:(n_generations + 1)) {
    cat("Seed:", seed, "| Generation:", gen, "\n")
    n <- nrow(individuals)
    reproduce <- runif(n) < 0.8
    reproducers <- individuals[reproduce, ]
    if (nrow(reproducers) == 0 && n > 0) {
      reproducers <- individuals[sample(seq_len(n), 1), ]
    }
    
    new_offspring <- list()
    min_dist <- 3 * circle_radius
    max_attempts <- 1
    
    for (i in seq_len(nrow(reproducers))) {
      parent <- reproducers[i, ]
      for (k in 1:5) {
        placed <- FALSE
        attempts <- 0
        while (!placed && attempts < max_attempts) {
          dispersal_sd <- circle_radius * 5
          dispersal_distance_x <- abs(rnorm(1, 0, dispersal_sd))
          dispersal_distance_y <- abs(rnorm(1, 0, dispersal_sd))
          x_new <- parent$x + dispersal_distance_x * sample(c(-1, 1), 1)
          y_new <- parent$y + dispersal_distance_y * sample(c(-1, 1), 1)
          if (x_new < 0 || x_new > space_size || y_new < 0 || y_new > space_size) break
          
          if (length(new_offspring) > 0) {
            xs <- vapply(new_offspring, function(o) o$x, numeric(1))
            ys <- vapply(new_offspring, function(o) o$y, numeric(1))
            dists <- sqrt((xs - x_new)^2 + (ys - y_new)^2)
            if (all(dists >= min_dist)) placed <- TRUE
          } else placed <- TRUE
          attempts <- attempts + 1
        }
        
        if (placed) {
          new_offspring[[length(new_offspring) + 1]] <- list(
            ID = NA, ParentID = parent$ID, Generation = gen,
            x = x_new, y = y_new, Sequence = parent$Sequence
          )
        }
      }
    }
    
    if (length(new_offspring) == 0) break
    
    offspring <- as.data.frame(do.call(rbind, lapply(new_offspring, as.data.frame)), stringsAsFactors = FALSE)
    offspring$ID <- next_ID:(next_ID + nrow(offspring) - 1)
    offspring[, c("ParentID", "Generation", "x", "y")] <- lapply(
      offspring[, c("ParentID", "Generation", "x", "y")], as.numeric)
    offspring$Sequence <- mutate_sequences(as.character(offspring$Sequence), mutation_rate)
    
    if (nrow(offspring) > K) offspring <- offspring[sample(nrow(offspring), K), ]
    
    offspring$sampling_year <- start_year + gen - 1
    offspring$Sequence_name <- paste0(offspring$ID, "_", offspring$sampling_year)
    
    all_individuals_list[[gen]] <- offspring
    individuals <- offspring
    next_ID <- next_ID + nrow(offspring)
    
    # Save files with seed
    write.fasta(as.list(offspring$Sequence), offspring$Sequence_name,
                sprintf("seed%d_all_individuals_sequences_gen%d.fasta", seed, gen), as.string = TRUE)
    write.table(offspring[, c("Sequence_name", "x", "y")],
                sprintf("seed%d_all_individuals_coordinates_gen%d.txt", seed, gen),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Plot
    p <- ggplot(offspring, aes(x = x, y = y)) +
      geom_rect(aes(xmin = 0, xmax = space_size, ymin = 0, ymax = space_size),
                fill = NA, color = "black", size = 1.2, inherit.aes = FALSE) +
      geom_point(color = "#e34a33", alpha = 0.7, size = 2) +
      coord_fixed(xlim = c(0, space_size), ylim = c(0, space_size)) +
      labs(title = paste("Seed", seed, "| Generation", gen, "| n =", nrow(offspring)),
           x = "X coordinate (normalized)",
           y = "Y coordinate (normalized)") +
      theme_minimal()
    print(p)
    cat("Finished generation", gen, "for seed", seed, "\n")
  }
}
