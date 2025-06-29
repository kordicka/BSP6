library(dplyr)
library(tidyverse)

# Set your current working directory containing the downloaded nodemap and interaction data from GeneGo MetaCore
# workingdir = '/Users/downloads/'  # MAC example
# workingdir = 'C:/Users/downloads/' # #Windows example
# setwd(workingdir) 
#

# set paths
data_dir <- "/home/nora/Semester6/BSP6/CODE/Network_Analysis/"
setwd(data_dir)

conditions = c("case_vs_controls")

for (condition in conditions) {
  
  #change and add "network" for nodes and interactions
  temp_nodes = paste0(condition, "_network_nodes")
  temp_interact = paste0(condition, "_network_interactions")
  
  # Reading network object and network interaction data from Excel files
  net_nodes = read.csv(paste0(data_dir, temp_nodes, ".csv"), 
                       skip = 2,     # Skip the first two metadata rows
                       sep = ",",    # Specify semicolon as separator
                       na.strings=c("","; ","NA"))
  net_interact = read.csv(paste0(data_dir, temp_interact, ".csv"), 
                          skip = 2,     # Skip the first two metadata rows
                          sep = ",",    # Specify semicolon as separator
                          na.strings=c("","; ","NA"))
  #We also transform each interaction "Unspecified"~"--?", "Activation"~"-->", "Inhibition" ~ "--|" so it can be used directly by Network pipeline
  net_nodes_processed = net_nodes %>% select(Network.Object.Name, Gene.Symbol) %>% drop_na()  %>% distinct() #We want  to remove duplicate and row containg NA
  write.table( net_nodes_processed , file =  paste0(temp_nodes, ".txt"), sep = "\t",row.names = FALSE, col.names = F, quote = F)
  
  #We also transform each interaction "Unspecified"~"--?", "Activation"~"-->", "Inhibition" ~ "--|" so it can be used directly by Network pipeline
  Interaction_processed= net_interact %>% select(Network.Object..FROM., Network.Object..TO., Effect) %>% drop_na() %>% filter(!(Effect=="Technical")) %>% mutate(interaction = case_match(Effect,"Unspecified"~"--?", "Activation"~"-->", "Inhibition" ~ "--|",.default = "ERROR")) %>% select(Network.Object..FROM., interaction, Network.Object..TO.)
  write.table(Interaction_processed , file =  paste0(temp_interact, ".txt"), sep = "\t",row.names = FALSE, col.names = F, quote = F)
  
  
}





