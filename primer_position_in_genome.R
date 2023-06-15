## Look for the primers of Mpox in some genomes ##
library(Biostrings)
library(tidyverse)

genome <- readDNAStringSet(filepath = "")
genome <- genome$`CP034122.1 Caulobacter vibrioides strain CB2A chromosome, complete genome`
genome <- as.character(genome)
primer_f <- "GGAAAGTGTAAAGACAACGAATACAG"
primer_r <- "CGATGTGGAAATTAACCTGTATC"


possible_attachments <- function (forward,reverse,genome,preferred_minimum_length,product_size){
  #this code assumes that the polymerase needs at least a primer of length 8 to have enough space to start polymerizing.
  #the user can change the preferred length and assume that the polymerase binds to a shorter fragment.
  trim_repetitions <- function(x){
    x <- x[x > 0]
    return(x)
  }
  master_list <- list()
    repetitions_f <- c()
    for (i in 1:nchar(forward)){
      substringed <- substr(forward,start=nchar(forward)-i+1,stop=(nchar(forward)))
      how_many <- lengths(regmatches(genome,
                                     gregexpr(pattern = substringed,text = genome)))
      repetitions_f[i] <- how_many
      if (how_many == 0){
        break
      }
    }
    repetitions_f <- trim_repetitions(repetitions_f)
    repetitions_r <- c()
    for (i in 1:nchar(reverse)){
      substringed <- substr(reverse,start=1,stop=i+1-1)
      how_many <- lengths(regmatches(genome,
                                     gregexpr(pattern = substringed,text = genome)))
      repetitions_r[i] <- how_many
      if (how_many == 0){
        break
      }
    }
    repetitions_r <- trim_repetitions(repetitions_r)
    if (missing(preferred_minimum_length) & length(repetitions_f) < 8 & length(repetitions_r) < 8){
      message("There are no matches with at least 8 bases. Specify a smaller annealing threshold")
    } else {
      if (missing(preferred_minimum_length)){
        threshold <- 8
      } else {
        threshold <- preferred_minimum_length
      }
    }
    repetitions_f <- c(threshold:length(repetitions_f))
    repetitions_r <- c(threshold:length(repetitions_r))
    combination_primers <- data.frame(matrix(ncol = 3))
    combination_primers_r <- data.frame(matrix(ncol = 3))
    #
    #
    for (i in 1:length(repetitions_f)){
      temporary_matrix <- matrix()
      string_interest <- substr(forward,start=nchar(forward)-repetitions_f[i]+1,stop=nchar(forward))
      h <- length(regmatches(genome,
                             gregexpr(string_interest,genome))[[1]])
      combination_primers <- rbind(combination_primers,
        data.frame(cbind(rep(string_interest,h),gregexpr(string_interest,genome)[[1]][1:h],rep("forward",h))))
    }
    #
    for (i in 1:length(repetitions_r)){
      temporary_matrix <- matrix()
      string_interest <- substr(reverse,start=1,stop=repetitions_r[i])
      h <- length(regmatches(genome,
                             gregexpr(string_interest,genome))[[1]])
      combination_primers_r <- rbind(combination_primers_r,
                                   data.frame(cbind(rep(string_interest,h),
                                                    gregexpr(string_interest,genome)[[1]][1:h],rep("reverse",h))))
    }
    #
    combination_primers <- combination_primers[-1,]
    combination_primers$X4 <- nchar(combination_primers$X1)+as.double(combination_primers$X2)
    combination_primers$X5 <- nchar(combination_primers$X1)
    filter_data_Frame <- data.frame(as_tibble(combination_primers) %>% 
                                      group_by(X4) %>% 
                                      summarise(max(X5)))
    colnames(filter_data_Frame)[2] <- "X5"
    a <- as_tibble(combination_primers) %>% 
      right_join(y = as_tibble(filter_data_Frame),by = c("X4","X5")) %>% 
      select(-c(X4,X5))
    #
    combination_primers_r <- combination_primers_r[-1,]
    combination_primers_r$X4 <- nchar(combination_primers_r$X1)
    filter_data_Frame_r <- data.frame(as_tibble(combination_primers_r) %>% 
                              group_by(X2) %>% 
                              summarise(max(X4)))
    colnames(filter_data_Frame_r)[2] <- "X4"
    b <- as_tibble(combination_primers_r) %>% 
      right_join(y = as_tibble(filter_data_Frame_r),by = c("X2","X4")) %>% 
      select(-X4)
    master <- data.frame(rbind(a,b))
    colnames(master) <- c("Sequence","Position","Orientation")
    master_list[[1]] <- master
    #
    #
    experimentals <- data.frame(matrix(ncol=5))
    col_names <- c("Sequence","Position","Orientation","product_size","V5")
    colnames(experimentals) <- col_names
    for (i in 1:nrow(master)){
      v <- abs(as.double(master[i,2])-as.double(master[,2]))
      experimental <- master
      experimental[,4] <- v
      if (sum((experimental$V4 < product_size),na.rm=T) > 1){
        experimental <- experimental[experimental$V4 < product_size,]
        colnames(experimental)[4] <- "product_size"
        experimental <- experimental[experimental$product_size > 0,]
        experimental[,5] <- paste0("with ",master[1,3]," ",master[1,1])
        experimental <- rbind(experimentals,experimental)} else {
          experimental <- data.frame(matrix(ncol=5))
          colnames(experimental)<-col_names
          experimental <- rbind(experimentals,experimental)
          next
      }
    }
    experimental <- experimental[!is.na(experimental[,5]),]
    if (nrow(experimental) == 0 ){
      master_list[[2]] <- "No products found"
    } else {
      master_list[[2]] <- experimental
    }
    return(master_list)
}

master <- possible_attachments(forward = primer_f,reverse = primer_r,genome = genome,product_size = 1500,preferred_minimum_length = 5)
as_tibble(master[[2]])


