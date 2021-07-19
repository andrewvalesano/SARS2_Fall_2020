

### Author: Andrew Valesano
### Purpose: Gather and analyze TreeTime data for ancestral state reconstruction with binary model (student vs. non-student).


# ============================= Modules and tree data ==========================

library(tidyverse)
library(ggtree)
library(ape)
library(wesanderson)
require(treeio)
require(phytools)
library(lubridate)
library(reshape2)
library(readxl)

# =============================== TreeTime output =================================

tree <- read.nexus("mugration_output/annotated_tree.nexus")

# TreeTime output: tip/node dates with confidence intervals
dates <- read.table("treetime_output/dates.tsv", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t") %>%
  rename(name = V1, date = V2, numeric_date = V3, lower_bound = V4, upper_bound = V5)

# Geographical regions used in mugration analysis and TreeTime inferences
regions <- c("Non-Student", "Student")
letters <- c("A", "B")
regions_df <- data.frame(letter = letters, region = regions)

read.csv("mugration_output/confidence.csv", stringsAsFactors = FALSE) %>%
  rename(name = X.name) %>%
  melt() %>% 
  rename(letter = variable, confidence = value) -> region_calls

region_calls <- left_join(region_calls, regions_df, by = "letter")
region_and_dates <- left_join(region_calls, dates, by = "name")

# Tip names and numbers
tips <- data.frame(tip_name = tree$tip.label)
tips <- mutate(tips, tip_num = 1:nrow(tips))

# Node names and numbers
nodes <- data.frame(node_name = tree$node.label)
nodes <- mutate(nodes, node_num = 1:nrow(nodes) + nrow(tips))

# Combine tips and nodes into a single dataframe
tree_data <- rbind(rename(tips, name = tip_name, num = tip_num) %>% mutate(type = "tip"), 
                   rename(nodes, name = node_name, num = node_num) %>% mutate(type = "node"))

# This dataframe has all of the treetime output (except for the tree itself).
treetime_data <- left_join(region_and_dates, tree_data, by = "name")

# All of the paths in the tree, in a list by tip number
paths <- nodepath(tree)

# Get tree paths to all tips
tip_paths <- data.frame()
for(t in tips$tip_num)
{
  path <- paths[[t]]
  path_df <- data.frame(step = path, tip_num = t)
  tip_paths <- rbind(tip_paths, path_df)
}

tip_paths <- mutate(tip_paths, step_type = ifelse(step == tip_num, "tip", "node"))
treetime_data <- mutate(treetime_data, step = num)
tip_paths_treetime <- left_join(tip_paths, treetime_data, by = "step") ### This is the working dataset.

# ================================= Find state transitions: Logic functions ============================

GetTransitionsToStudent <- function(df, confidence_level)
{
  df_transition <- filter(df, region == "Student" & confidence >= confidence_level)
  
  if(nrow(df_transition) == 0)
  {
    df_none <- mutate(df, region = "NONE") # Send this back so we aren't returning NULL dataframes, can filter later
    return(df_none)
  }
  else
  {
    transition_node <- df_transition[1,]
    return(transition_node)
  }
}

# For paths containing Student nodes, do any go back to Non-Student?
GetTransitionsFromStudents <- function(df, confidence_level)
{ 
  # If we passed a path that doesn't go through Student, we need to leave
  df_student <- filter(df, region == "Student" & confidence >= confidence_level)
  if(nrow(df_student) == 0)
  {
    df_none <- mutate(df, region = "NONE") # Send this back so we aren't returning NULL dataframes, can filter later
    return(df_none)
  }
  
  # Assign numbers to make it easy
  local_stepnums <- rep(seq(1, length(unique(df$step))), each=2)
  df <- mutate(df, local_stepnum = local_stepnums)
  
  # Find the first time we find "Student" confidently
  first_student_node_localstep <- 10000000000
  for(l in 1:max(df$local_stepnum))
  {
    df_stepnum <- filter(df, local_stepnum == l)
    df_stepnum_student <- filter(df_stepnum, region == "Student")
    df_stepnum_student_confidence <- unique(df_stepnum_student$confidence)
    if(df_stepnum_student_confidence > confidence_level & l < first_student_node_localstep)
    {
      first_student_node_localstep <- l
    }
  }

  # We reached a tip, not a node. Time to leave
  if(first_student_node_localstep == max(df$local_stepnum))
  {
    df_none <- mutate(df, region = "NONE") 
    return(df_none)
  }
  
  # Starting from the first "Student" node, find the first confident Non-Student
  df_afterstudent <- filter(df, local_stepnum > first_student_node_localstep)
  df_nonstudent <- filter(df_afterstudent, region == "Non-Student" & confidence >= confidence_level)
  transition_node <- filter(df_nonstudent, local_stepnum == min(local_stepnum))
  return(transition_node)
}

# ================================= Find state transitions ============================

confidence_levels <- c(0.9, 0.99)
nodes_to_students_all <- data.frame()
nodes_from_students_all <- data.frame()

for(c in confidence_levels)
{
  ### Get transitions from non-student to student
  tip_paths_treetime %>%
    group_by(tip_num) %>%
    do(GetTransitionsToStudent(., c)) %>%
    ungroup() %>%
    filter(!region == "NONE") -> nodes_to_students
  
  nodes_to_students_distinct <- mutate(nodes_to_students, confidence_level = c) %>%
    select(-tip_num) %>%
    distinct() %>%
    mutate(intro_type = ifelse(step_type == "node", "Non-Singleton", "Singleton"))
  
  nodes_to_students_all <- rbind(nodes_to_students_all, nodes_to_students_distinct) # Store for later
  
  ### Get transitions from students back to non-students
  student_nodes <- filter(nodes_to_students_distinct, step_type == "node")
  tip_paths_treetime_student_steps <- filter(tip_paths_treetime, step %in% student_nodes$step)
  tip_paths_treetime_student_paths <- filter(tip_paths_treetime, tip_num %in% tip_paths_treetime_student_steps$tip_num) # Paths containing student nodes
  
  tip_paths_treetime_student_paths %>%
    group_by(tip_num) %>%
    do(GetTransitionsFromStudents(., c)) %>%
    ungroup() %>%
    filter(!region == "NONE") -> nodes_from_students
  
  nodes_from_students_distinct <- mutate(nodes_from_students, confidence_level = c) %>%
    select(-tip_num) %>%
    distinct() %>%
    mutate(intro_type = ifelse(step_type == "node", "Non-Singleton", "Singleton"))
  
  nodes_from_students_all <- rbind(nodes_from_students_all, nodes_from_students_distinct)
}


# ================================== For these nodes, get descendant sizes and plot =============================

# Data for all confidence levels: nodes_to_students_all, nodes_from_students_all.

# Select a confidence level cutoff (one of the above)
cl <- 0.99

nodes_to_students <- filter(nodes_to_students_all, confidence_level == cl)
nodes_from_students <- filter(nodes_from_students_all, confidence_level == cl)

# Raw number of non-singleton transitions
nrow(filter(nodes_to_students, step_type == "node"))
nrow(filter(nodes_from_students, step_type == "node"))

GetDescendants <- function(node_data, descendant_type) 
{
  intro_type <- unique(node_data$intro_type)
  node <- unique(node_data$name)
  if(intro_type == "Singleton")
  {
    return_data <- data.frame(name = node, num_descendants = 1, num_descendants_other = 0)
    return(return_data)
  }
  else if(intro_type == "Non-Singleton")
  {
    tree.node <- tree_subset(tree, node = node, group_name = "subset_group", levels_back = 0)
    node_tips <- filter(treetime_data, name %in% tree.node$tip.label & region %in% descendant_type)
    tips <- unique(node_tips$name)
    tips_group <- tips[tips %in% filter(metadata, student_status == descendant_type)$strain]
    num_descendants <- length(tips_group)
    tips_group_other <- tips[tips %in% filter(metadata, student_status != descendant_type)$strain]
    num_descendants_other <- length(tips_group_other)
    return_data <- data.frame(name = node, num_descendants = num_descendants, num_descendants_other = num_descendants_other)
    return(return_data)
  }
  else
  {
    print("Houston we have an error!")
  }
}

nodes_to_students %>%
  group_by(name) %>%
  do(GetDescendants(., c("Student"))) %>%
  mutate(type = "To students") -> transitions_to_students

nodes_from_students %>%
  group_by(name) %>%
  do(GetDescendants(., c("Non-Student"))) %>%
  mutate(type = "From students") -> transitions_from_students

### Plot spillover size by introduction size ###

transitions_to_students_filtered <- filter(transitions_to_students, num_descendants > 1)

descendants.by.size.plot <- ggplot(transitions_to_students_filtered, aes(x = num_descendants, y = num_descendants_other)) +
  geom_point(size = 3, shape = 21, fill = "#75988d") +
  theme_bw(base_size = 15) +
  xlim(c(0, 120)) +
  ylim(c(0, 30)) +
  xlab("Size of Student Cluster") +
  ylab("Number of Non-Student Descendants") +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5)
  

