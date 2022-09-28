if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("magrittr", quietly = TRUE))
  install.packages("magrittr")

if (!require("readr", quietly = TRUE))
  install.packages("readr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

# Install the Human package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

expression_df <- readr::read_tsv("GSE183516_Normalized_counts_AT2_pAT2_samples.txt")

# Tuck away the Gene ID column as row names
expression_df = tibble::column_to_rownames(expression_df,var="...1")

expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)

# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Entrez") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Entrez)

multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "entrez_id_count") %>%
  # Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(entrez_id_count))

# Let's look at the first 6 rows of our `multi_mapped` object
head(multi_mapped)

expression_df_data <- expression_df
for( i in 1:nrow(expression_df_data)){
  expression_df_data[i, 1] <- mapped_df[i, 2] 
}

df_data_var <- expression_df_data
for( i in 1:nrow(df_data_var)){
  df_data_var[i, 2] <- var(as.vector(t(expression_df_data[i,2:ncol(df_data_var)])))
  #print(as.vector(t(expression_df_data[i,2:ncol(df_data_var)])))
}
geneVariation <- df_data_var[,1:2]
#print(as.vector(t(expression_df_data[1,2:ncol(df_data_var)])))

library(ggplot2)
p <- ggplot(geneVariation, aes(x = Gene, y = MUG207A37)) + 
  geom_point()
view(p)


# Log base 10 scale + log ticks (on left and bottom side)
p + scale_y_continuous(trans = 'log2')
view(p)
