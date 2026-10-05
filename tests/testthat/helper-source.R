source_files <- list.files(file.path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
for (source_file in source_files) source(source_file)

