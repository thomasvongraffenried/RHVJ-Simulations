library(seqinr)
library(ape)
library(ggplot2)

# -----------------------------
# USER PARAMETERS
# -----------------------------

base_dir <- "C:/Users/tv19s091/OneDrive - Universitaet Bern/Dokumente/Master/Master_Thesis/Data/HepJ/Simulation_23.12.25/LDD"

seeds <- 1:27
first_gen <- 1
last_gen <- 300
n_subsample <- 25

area_defs <- data.frame(
  Area = c("A", "B"),
  xmin = c(2/3, 0),
  xmax = c(1, 1/3),
  ymin = c(2/3, 0),
  ymax = c(1, 1/3),
  fill = c("red", "blue"),
  stringsAsFactors = FALSE
)

set.seed(1)

# -----------------------------
# LOOP OVER SEEDS
# -----------------------------

for (seed in seeds) {
  
  cat("\n========================================\n")
  cat("Searching for seed", seed, "\n")
  cat("========================================\n")
  
  found_gen <- NA
  
  for (gen in first_gen:last_gen) {
    
    coord_file <- file.path(
      base_dir,
      sprintf("seed%d_all_individuals_coordinates_gen%d.txt", seed, gen)
    )
    seq_file <- file.path(
      base_dir,
      sprintf("seed%d_all_individuals_sequences_gen%d.fasta", seed, gen)
    )
    
    if (!file.exists(coord_file) || !file.exists(seq_file)) next
    
    coords <- read.table(coord_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    seqs <- read.fasta(seq_file, as.string = FALSE, forceDNAtolower = FALSE)
    
    coords$Sequence_name <- trimws(coords$Sequence_name)
    names(seqs) <- trimws(names(seqs))
    idx <- match(coords$Sequence_name, names(seqs))
    if (any(is.na(idx))) next
    
    enough <- TRUE
    for (i in seq_len(nrow(area_defs))) {
      area <- area_defs[i, ]
      n_in_area <- sum(
        coords$x >= area$xmin & coords$x <= area$xmax &
          coords$y >= area$ymin & coords$y <= area$ymax
      )
      if (n_in_area < n_subsample) {
        enough <- FALSE
        break
      }
    }
    
    if (enough) {
      found_gen <- gen
      break
    }
  }
  
  if (is.na(found_gen)) {
    cat("❌ No suitable generation found for seed", seed, "\n")
    next
  }
  
  cat("✅ Seed", seed, ": suitable generation found =", found_gen, "\n")
  
  # -----------------------------
  # ASK USER
  # -----------------------------
  
  answer <- readline(
    sprintf("Do you want to subsample seed %d, generation %d? (yes/no): ",
            seed, found_gen)
  )
  
  if (tolower(answer) != "yes") {
    cat("➡️  Skipping seed", seed, "\n")
    next
  }
  
  # -----------------------------
  # LOAD DATA
  # -----------------------------
  
  coord_file <- file.path(
    base_dir,
    sprintf("seed%d_all_individuals_coordinates_gen%d.txt", seed, found_gen)
  )
  seq_file <- file.path(
    base_dir,
    sprintf("seed%d_all_individuals_sequences_gen%d.fasta", seed, found_gen)
  )
  
  coords <- read.table(coord_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  seqs <- read.fasta(seq_file, as.string = FALSE, forceDNAtolower = FALSE)
  
  coords$Sequence_name <- trimws(coords$Sequence_name)
  names(seqs) <- trimws(names(seqs))
  idx <- match(coords$Sequence_name, names(seqs))
  seqs <- seqs[idx]
  
  # -----------------------------
  # SUBSAMPLING
  # -----------------------------
  
  sub_idx <- list()
  already_sampled <- integer(0)
  
  for (i in seq_len(nrow(area_defs))) {
    area <- area_defs[i, ]
    idx_area <- which(
      coords$x >= area$xmin & coords$x <= area$xmax &
        coords$y >= area$ymin & coords$y <= area$ymax
    )
    available <- setdiff(idx_area, already_sampled)
    chosen <- sample(available, n_subsample)
    
    sub_idx[[area$Area]] <- chosen
    already_sampled <- c(already_sampled, chosen)
  }
  
  # -----------------------------
  # WRITE FILES
  # -----------------------------
  
  coords_list <- list()
  seqs_list <- list()
  
  for (area_name in names(sub_idx)) {
    
    coords_area <- coords[sub_idx[[area_name]], ]
    coords_area$Sequence_name <- paste0(coords_area$Sequence_name, "_", area_name)
    
    seqs_area <- seqs[sub_idx[[area_name]]]
    names(seqs_area) <- coords_area$Sequence_name
    
    coords_list[[area_name]] <- coords_area
    seqs_list[[area_name]] <- lapply(seqs_area, as.character)
    
    write.fasta(
      seqs_list[[area_name]],
      names(seqs_list[[area_name]]),
      file.path(
        base_dir,
        sprintf("seed%d_subsample_%s_gen%d.fasta", seed, area_name, found_gen)
      )
    )
    
    write.table(
      coords_area,
      file.path(
        base_dir,
        sprintf("seed%d_subsample_%s_gen%d_coordinates.txt", seed, area_name, found_gen)
      ),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  
  coordsAB <- do.call(rbind, coords_list)
  seqsAB <- do.call(c, seqs_list)
  names(seqsAB) <- coordsAB$Sequence_name
  
  write.fasta(
    seqsAB, names(seqsAB),
    file.path(base_dir,
              sprintf("seed%d_subsample_AB_gen%d.fasta", seed, found_gen))
  )
  
  ape::write.nexus.data(
    seqsAB,
    file = file.path(base_dir,
                     sprintf("seed%d_subsample_AB_gen%d.nexus", seed, found_gen)),
    format = "dna"
  )
  
  # -----------------------------
  # PLOT (SAME STYLE AS ORIGINAL)
  # -----------------------------
  
  coords$Group <- "All"
  for (area_name in names(sub_idx)) {
    coords$Group[sub_idx[[area_name]]] <- paste("Area", area_name)
  }
  
  p <- ggplot(coords, aes(x = x, y = y)) +
    geom_point(color = "gray70", size = 2) +
    geom_point(data = coords[coords$Group == "Area A", ],
               color = "red", size = 4) +
    geom_point(data = coords[coords$Group == "Area B", ],
               color = "blue", size = 4) +
    geom_rect(
      data = area_defs,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Area),
      color = "black", alpha = 0.2, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = setNames(area_defs$fill, area_defs$Area)) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = sprintf("Seed %d — Subsampling Generation %d", seed, found_gen),
      subtitle = "Rectangles show sampling areas; colored points are sampled individuals",
      x = "X coordinate (normalized)",
      y = "Y coordinate (normalized)",
      fill = "Area"
    ) +
    theme_minimal()
  
  print(p)
  
  ggsave(
    filename = file.path(
      base_dir,
      sprintf("seed%d_subsampling_gen%d_plot.png", seed, found_gen)
    ),
    plot = p,
    width = 7,
    height = 7,
    dpi = 300
  )
  
  cat("✅ Finished seed", seed, "\n")
}
