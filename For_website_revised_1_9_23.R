
# app.R


column(12, offset=3,align="center", textOutput('Purpose_text'), # 1/5/23
     tags$head(tags$style("#Purpose_text{color: black;
                                 font-size: 24px;
                                 font-style: serif;
                                 }"
                         )
     )
     ),

column(12, offset=3,align="left", textOutput('Steps_text'), # 1/5/23
     tags$head(tags$style("#Steps_text{color: blue;
                                font-size: 18px;
                                font-style: italic;
                                 }"
                         )
     )
     ),

##1/5/23
column(12, offset=3,align="left", h4("Step1.")),
column(12, offset=3,align="left", h5("A). Input the gene model ID in the Gene Name box, then click “(1) Check Gene ID” button to fill the Boxes for Reference and Chromosome.")),
column(12, offset=3,align="left", h5("B). To submit a transcript sequence, manually select the Reference and Chromosome, select “fasta_seq” in the “CDS” box, then past the sequence in the “Your fasta sequence” box.")),
column(12, offset=3,align="left", h5("C). Modify the value for Upstream and Downstream box to set the search boundaries. Select or unselect the Genomes to be included in the analysis. The click “(2) Submit) button to start the process; the alignment with all selected genomes will be plotted in the top right panel. Note: to upload your own chromosome into the analysis, click “Browse” button.")),

column(12, offset=3,align="left", h4("Step2.")),
column(12, offset=3,align="left", h5("Click the “(3) Tree” button to plot the phylogenetic tree clustering genomes based on shared indels. Click the tree plot (top right) to determine the number of haplotypes.")),

column(12, offset=3,align="left", h4("Step3.")),
column(12, offset=3,align="left", h5("Click the “(4) Plot selected haplotypes” button to plot the alignment among genomes representing each haplotype. To trim the unwanted regions from the final haplotype presentation, single click on the third panel to set the left boundary, then double click to set the right boundary.")),

column(12, offset=3,align="left", h4("Step4.")),
column(12, offset=3,align="left", h5("If coordinates of both boundaries are set, the  “(5) Trim” button will be clickable.")),
##1/5/23


column(12, offset=3,align="left", textOutput('Example_text'), # 1/5/23
     tags$head(tags$style("#Example_text{color: blue;
                                 font-size: 18px;
                                 font-style: italic;
                                 }"
                         )
     )
     ),




column(12, offset=3,align="left", textOutput('Algorithm_text'), # 1/5/23
     tags$head(tags$style("#Algorithm_text{color: blue;
                                 font-size: 18px;
                                 font-style: italic;
                                 }"
                         )
     )
     ),
column(12, offset=3,align="center", tags$img(width="600", height="400", src="https://bridgecereal.scinet.usda.gov/BRIDGEcereal_Default_Parameters.png")), # 1/5/23




column(12, offset=3,align="left", textOutput('Demo_text'), # 1/5/23
     tags$head(tags$style("#Demo_text{color: blue;
                                 font-size: 18px;
                                 font-style: italic;
                                 }"
                         )
     )
     ),




column(12, offset=3,align="center", h5("Acknowledgements: We thank the USDA-ARS SCINet for computing resource and the collaboration of the USDA-ARS-Partnerships for Data Innovations (PDI, https://pdi.scinet.usda.gov/), which provided data stewardship solutions to enable secure data management, storage and sharing. Contact: xianran.li@usda.gov")),# 1/5/23
column(5,offset=6, align="center", tags$img(width="120", height="40", src="https://bridgecereal.scinet.usda.gov/USDA_PDI_Logo.jpg")), # 1/5/23


### function part

output$Purpose_text <- renderText({                                 # 1/5/23
paste("BRIDGEcereal: An interactive webapp to discover indel variation in cereal pan-genomes",sep="")
  })

output$Steps_text <- renderText({                                 # 1/5/23
paste("BRIDGEcereal app brief instructions in 4 steps:",sep="")
  })

output$Example_text <- renderText({                                 # 1/5/23
paste("Figure1: Instructions in 4 steps (highlighted by yellow color) on one wheat (gene) example:",sep="")
  })

output$Algorithm_text <- renderText({                                 # 1/5/23
paste("Figure2: Algorithm used to determine query gene's orthologs (after clicking the Submit button):",sep="")
  })

output$Demo_text <- renderText({                                 # 1/5/23
paste("BRIDGEcereal app short viedo instructions:",sep="")
  })



### sub.R, For plot4 only

#output_flag = 0 # 1/9/23
output_flag = 1 # 1/9/23


output_flag = 0 # 1/9/23
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, output_flag,Gene,Ref_genome,haplotypes,repeats,strand_direction) # 1/9/23


#BlastSynWorking_2 <- BlastSynWorking[which(BlastSynWorking$Genome!=''),] # 1/9/23
BlastSynWorking_2 <- Filtered_HaplotypeSyn_Plotted # 1/9/23

