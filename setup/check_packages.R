pkgs <- c('doSNOW','doParallel','doMPI','tidyverse','tictoc','fields','ggplot2',
          'sp','lattice','plotly','geometry','viridis','gstat','remotes',
          'kableExtra','splancs','rgdal','INLAOutputs','inlamesh3d',
          'inlatools','INLAutils','ggregplot','INLA')
installed <- rownames(installed.packages())
missing <- pkgs[!pkgs %in% installed]
cat("Installed:\n")
cat(paste(pkgs[pkgs %in% installed], collapse="\n"), "\n\n")
cat("Missing:\n")
cat(paste(missing, collapse="\n"), "\n")
