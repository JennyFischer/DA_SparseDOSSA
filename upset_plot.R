library(UpSetR)
listInput <- list(
                  Aldex2_Wilcoxon_Relative_Abundance_CLR = c(96,99,155),
                  Aldex2_TTest_Relative_Abundance_CLR = c(),
                 TTest_Relative_Abundance_CSS= c(19,21,33,42,45,78),
                 TTest_Relative_Abundance_TSS= c(19,21,33,42,45,78),
                 TTest_Relative_Abundance_CLR= c(19,21,33,42,45,78),
                 TTest_Relative_Abundance_RA= c(19,21,33,42,45,78),
                 Wicoxon_Relative_Abundance_CSS= c(19,21,33,42,45,78),
                 Wicoxon_Relative_Abundance_TSS= c(19,21,33,42,45,78),
                Wicoxon_Relative_Abundance_CLR= c(19,21,33,42,45,78),

                Wicoxon_Relative_Abundance_RA= c(19,21,33,42,45,78),
 
                Wilcoxon_Absolut_Abundance = c(7,19,21,33,35,42,43,45,48,79))
listInput          

upset(fromList(listInput),nsets= 11)



listInput_2set <- list(
                  TTest_Relative_Abundance_CSS= c(129,144,151,218,234),
                  Aldex2_TTest_Relative_Abundance_CLR = c(),
                  Aldex2_Wilcoxon_Relative_Abundance_CLR = c(3,109, 139, 238),
                  TTest_Relative_Abundance_TSS= c(129,144,151,218,234),
                  TTest_Relative_Abundance_CLR= c(129,144,151,218,234),
                  TTest_Relative_Abundance_RA=c(129,144,151,218,234),
                  Wicoxon_Relative_Abundance_CSS= c(129,144,151,218,234),
                  Wicoxon_Relative_Abundance_TSS= c(129,144,151,218,234),
                  Wicoxon_Relative_Abundance_CLR= c(129,144,151,218,234),
                  Wicoxon_Relative_Abundance_RA= c(129,144,151,218,234),
                  Wilcoxon_Absolut_Abundance= c(129,144,151,218,234))

upset(fromList(listInput_2set),nsets= 11)

