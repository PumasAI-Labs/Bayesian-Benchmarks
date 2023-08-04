lowTriMat <- function(x) {
  diag(x)[upper.tri(diag(x), diag = T)]
}

get_theta <- function(run, chain){
  
  read_json(str_c(root_folder, "/data/inits/inits_", run, "_", chain, 
                  ".json")) %>% 
    spread_all() %>% 
    as_tibble() %>% 
    select(starts_with("TV")) %>% 
    mutate(across(everything(), log),
           across(everything(), \(x) round(x, 3))) %>% 
    as.numeric() %>% 
    stringi::stri_paste(collapse = " ")
  
}

get_omega <- function(run, chain){
  
  omegaMat <- read_json(str_c(root_folder, "/data/inits/inits_", run, "_", 
                              chain, ".json")) %>% 
    enter_object(omega) %>%
    gather_array() %>%
    # spread_all() %>%
    pull(..JSON) %>%
    unlist() %>% 
    lowTriMat() 
  
  omegaMat[omegaMat == 0] <- 0.01
  
  omegaMat %>% 
    stringi::stri_paste(collapse = " ")
  
}

get_sigma <- function(run, chain){
  
  read_json(str_c(root_folder, "/data/inits/inits_", run, "_", chain, 
                  ".json")) %>% 
    spread_all() %>% 
    as_tibble() %>% 
    select(sigma) %>% 
    as.character()
  
}

write_text_line <- function(string = c("ITERFILE", "SDFILE", "PAFILE"), run,
                            chain, modelText){
  
  string <- match.arg(string)
  
  type <- case_when(string == "ITERFILE" ~ str_c("-iter-", "run", run, "-"),
                    string == "SDFILE" ~ str_c("-sdtab-", "run", run, "-"),
                    string == "PAFILE" ~ str_c("-patab-", "run", run, "-"),
                    .default = NA_character_)
  
  modelTextLine <- grep(string, modelText)
  
  modelText[modelTextLine] <- gsub(string,  
                                   paste0(modelName, type, chain, ".tab"), 
                                   modelText[modelTextLine])
  
  return(modelText)
  
}

create_model_file <- function(run, chain){
  
  modelText <- readLines(paste0(dir_name, "/", modelName, "-template.mod"))
  
  modelTextLine <- grep("THETAUPDATE", modelText)
  modelText[modelTextLine] <- get_theta(run, chain)
  
  modelTextLine <- grep("OMEGAUPDATE", modelText)
  modelText[modelTextLine] <- get_omega(run, chain)
  
  modelTextLine <- grep("SIGMAUPDATE", modelText)
  modelText[modelTextLine] <- get_sigma(run, chain)
  
  modelText <- write_text_line("ITERFILE", run, chain, modelText)
  modelText <- write_text_line("SDFILE", run, chain, modelText)
  modelText <- write_text_line("PAFILE", run, chain, modelText)
  
  write(modelText, paste0(base_name, "-run", run, "-", chain,".mod"))
  
}