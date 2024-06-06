library(tidyverse)
random <- read.csv('~/Desktop/GenotypesTreatments.csv')

#
sampling <- random %>%
  select(CoralNum, Genotype)%>%
  group_by(Genotype) %>%
  mutate(ID = sample(1:4, 4, replace = FALSE)) %>%
  ungroup()%>%
  mutate(Treatment = case_when(ID == 1 ~ "PRE",
                               ID == 2 ~ "COTS",
                               ID == 3 ~ "BLEACH",
                               ID == 4 ~ "CONTROL"))
sampling

write.csv(x = sampling, file = '~/Desktop/corals_random_assignments.csv')

