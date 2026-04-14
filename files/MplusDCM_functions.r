#Functions for analysis- (by Jonathan Templin)-----------------------------------------------------------------------------------
dec2bin = function(decimal_number, nattributes, basevector){
  dec = decimal_number
  profile = matrix(NA, nrow = 1, ncol = nattributes)
  for (i in nattributes:1){
    profile[1,i] =  dec%%basevector[1,i]
    dec = (dec-dec%%basevector[1,i])/basevector[1,i]
  }
  return(profile)
}

marginalize_attribute_probabilities = function(profile_estimates, natts){
  #create attribute profile matrix
  att_profile = NULL
  for (profile in 0:(2^natts-1)){
    att_profile = rbind(att_profile, dec2bin(decimal_number = profile, nattributes = 3, basevector = matrix(rep(2,3), nrow = 1)))
  }
  
  #marginalize attributes via matrix product
  marginal_estimates = as.matrix(profile_estimates) %*% att_profile
  
  marginal_estimates = data.frame(marginal_estimates)
  colnames(marginal_estimates) = paste("Attribute", 1:natts, sep="")
  return(marginal_estimates)
}

calculate_attribute_reliability = function(marginal_attribute_estimates, nreps){
  
  if (tryCatch(require("polycor")) == FALSE){
    install.packages("polycor")
    library(polycor)
  }
  
  #convert to matrix:
  marginal_attribute_estimates = as.matrix(marginal_attribute_estimates)
  
  #draw probabilities with replication:
  sampled_examinees = as.matrix(sample(x = marginal_attribute_estimates, size = nreps, replace = TRUE))
  
  #draw data, twice
  simdata = cbind(rbinom(n=nreps, size=1, prob=sampled_examinees), rbinom(n=nreps, size=1, prob=sampled_examinees))
  
  #calculate tetrachoric correlation between two draws
  reliability = polychor(x = simdata[,1], y = simdata[,2])
  
  return(reliability)
}


Mplus_to_LCDM = function(qmatrix, MplusAnalysis){
  #get attribute names
  attribute_names = colnames(qmatrix)[2:dim(qmatrix)[2]]
  
  #get item names
  item_names = qmatrix$Item
  
  #extract LCDM parameters:
  LCDM_item_parameters = MplusAnalysis$parameters$unstandardized[which(MplusAnalysis$parameters$unstandardized$paramHeader == "New.Additional.Parameters"),]
  
  #remove any parameters with "g" in label (structural model):
  LCDM_item_parameters = LCDM_item_parameters[grep(pattern = "L", x = LCDM_item_parameters$param),]
  
  #add columns to parameters
  LCDM_item_parameters$mplus_name = LCDM_item_parameters$param
  
  #get item number
  i=2
  for (i in 1:dim(LCDM_item_parameters)[1]){
    item_split = strsplit(x = LCDM_item_parameters$mplus_name[i], split = "_")
    
    item_num_text = item_split[[1]][1]
    
    item_num = as.numeric(substr(x = item_num_text, start = 2, stop = nchar(item_num_text)))
    
    item_name = item_names[item_num]
    
    LCDM_item_parameters$item[i] = item_name
    
    effect_text = item_num_text = item_split[[1]][2]
    effect_level = as.numeric(substr(x = effect_text, start = 1, stop=1))
    LCDM_item_parameters$effect_level[i] = effect_level
    
    if (effect_level == 0){
      LCDM_item_parameters$effect[i] = "Intercept"
    } else if (effect_level == 1){
      LCDM_item_parameters$effect[i] = "Main Effect"
      effect_att = as.numeric(substr(x = effect_text, start = 2, stop=2))
      effect_att_name = attribute_names[effect_att]
      
      LCDM_item_parameters$effect_att[i] = effect_att_name
      
    } else if (effect_level >= 2){
      LCDM_item_parameters$effect[i] = "Interaction"
      
      effect_att_name = NULL
      for (j in 1:effect_level){
        effect_att = as.numeric(substr(x = effect_text, start = (j+1), stop=(j+1)))
        if (j==1){
          effect_att_name = attribute_names[effect_att]
        } else {
          effect_att_name = paste(effect_att_name, attribute_names[effect_att], sep="*")
        }
      }
      LCDM_item_parameters$effect_att[i] = effect_att_name
    }
    
  }
  
  LCDM_item_parameters = LCDM_item_parameters[c("mplus_name", "item", "effect_level", "effect", "effect_att", "est", "se", "est_se", "pval")]
  
  #next: extract attribute probabilities from savedata file
  
  #get number of attributes
  natts = length(attribute_names)
  npattern = 2^natts
  
  #create attribute profile matrix
  att_profile = NULL
  att_profile_names = NULL
  profile = 0
  for (profile in 0:(2^natts-1)){
    att_profile = rbind(att_profile, dec2bin(decimal_number = profile, nattributes = natts, basevector = matrix(rep(2,natts), nrow = 1)))
    temp_name = "P"
    for (att in 1:natts){
      temp_name = paste(temp_name, att_profile[(profile+1),att], sep="")
    }
    att_profile_names = c(att_profile_names, temp_name)
  }
  
  #create Mplus names
  MplusProbNames = paste("CPROB", 1:npattern, sep="")
  MplusClassName = "C"
  
  EAP_profile_estimates = MplusAnalysis$savedata[MplusProbNames]
  colnames(EAP_profile_estimates) = att_profile_names
  
  MAP_profile_estimate = MplusAnalysis$savedata[MplusClassName]
  MAP_profile_estimate$Profile = att_profile_names[MAP_profile_estimate$C]
  
  #create marginal attribute profile estimates 
  EAP_marginal_estimates = data.frame(as.matrix(EAP_profile_estimates) %*% att_profile)
  colnames(EAP_marginal_estimates) = attribute_names
  
  MAP_marginal_estimates = round(EAP_marginal_estimates, digits = 0)
  
  #create attribute reliabilities
  attribute_information = data.frame(attribute = attribute_names, stringsAsFactors = FALSE)
  attribute_information$reliability = NA
  
  for (i in 1:natts){
    attribute_information$reliability[i] = calculate_attribute_reliability(marginal_attribute_estimates = EAP_marginal_estimates[,i], nreps = 100000)  
  }
  
  return(list(LCDM_item_parameters = LCDM_item_parameters, attribute_reliability = attribute_information, EAP_marginal_estimates = EAP_marginal_estimates, 
              EAP_profile_estimates = EAP_profile_estimates, MAP_marginal_estimates = MAP_marginal_estimates, MAP_profile_estimates = MAP_profile_estimate))
}




 ############################################################################## 
 #              Latent Class Attribute Profile Generator    by Rupp and Wilhelm                  #
 ##############################################################################

 generateLatentClasses <- function(qmatrix){
 
      num_att = length(qmatrix[1,])
      max_att = 2^num_att
      latent_classes = matrix (data=NA, max_att, num_att)
      m <- max_att

      for (a in 1:num_att) {
          m = m/2     # Number of repititions of entries 0 or 1 in one cycle
          anf <- 1

          while (anf < max_att) {
            latent_classes[anf : (anf + m -1),a] <- 0
            anf <- anf + m
            latent_classes[anf : (anf + m -1),a] <- 1
            anf <- anf + m
          }
      }
      rownames(latent_classes) = paste("c", 1:max_att, sep= "")
      latent_classes
      
 } #### end of function generateLatentClasses
 
 

###############################################################################
#                           INPUT FILE GENERATOR                              #
###############################################################################

# Generates an input file for LCDM model for MPLUS
# incl. the required MODEL command and MODEL constrains
# Main effects as well as 1- till 5-way-interactions are implemented

# Input: Data set (only relevant items!) + Q-matrix + maximum level of interaction +
#        Rule for the LCDM (FULL, CRUM, NIDO)   
# Output: Input file for Mplus

# Note that in the produced input file lines can be longer than 80 characters.
# In this case MPlus will show a warning or an error massage. Shorten the lines by 
# hand before starting the model estimation in Mplus in this case.

################################################################################

generate.mplus.input <- function(filename_ds, qmatrix, input_file, rule="FULL", max_level_interaction=3,
                                 orderConstrains=T, PreAddVars="", PostAddVars="", miss=F, miss_char=7 ){

    # library(xlsReadWrite)
    
    # Number of attributes and items in the Q-matrix
    anzahl_att <- length(qmatrix[1,])
    anzahl_items <- length(qmatrix[,1])
        
    # If one of the names is 1 rename all items to x1
    # First try to deal with item names starting with a number
    if (1%in%row.names(qmatrix)) {
        row.names(qmatrix) <- paste("x", row.names(qmatrix), sep="")
    }
    
    # Adoption of the number of inteactions depending on the selected RULE
    if ((rule=="CRUM") || (rule=="NIDO") ) { max_level_interaction <- 1 }
       
    latent_classes <- generateLatentClasses(qmatrix)
    max_att = 2^anzahl_att
    
    #--------------------------------------------------------------------------------------------------
    #------- Create LCDM kernels and order restrictions for all Items ---------------------------------
    
    # generate matrix with kernel names and constrains - only for testing purposes
    kernel_namen = matrix(data=NA, anzahl_items, max_att)
    colnames(kernel_namen) <- rownames(latent_classes)
    
    kernel_item_constrains = matrix(data=NA, anzahl_items, max_att)
    colnames(kernel_item_constrains) <- rownames(latent_classes)
    
    i <- 1        # i - item
    a <- 1        # a, j, k, l - attributs
    lc <- 1       # lc -latent class
    
    # initialisation - contain later all kernals and order restrictions
    all_kernel_constrains = ""
    all_order_constrains = ""
    kernel_item_without_names =""
    NEW_params_Mplus = ""      # includes later all parameters for the MODEL CONSTRAIN NEW part 
    
    ####### loop over all items
    for (i in 1:anzahl_items) {
    
        m <- 0                 # iterator for the names of the kernals for the current item
        anz_att_items <- sum (qmatrix[i,])  
    
        ########## loop over all latent classes                
        for (lc in 1:max_att) {
          
          # basic parameter            
          kernel_item_constrain <- paramname <- paste("l", i, "_0", sep="")      
          if (paramname%in%NEW_params_Mplus == FALSE) 
               NEW_params_Mplus = c(NEW_params_Mplus, paramname)
              
          ################### Main effects for attributes
          if (anz_att_items > 0) {      
            for (a in 1:anzahl_att) {      
              if ((qmatrix[i,a]==1) & (latent_classes[lc,a]==1)) {
    
                # add the main effects depending on the selected rule            
                if (rule!="NIDO") {
                      paramname <- paste ("l", i, "_1", a ,sep="")
                }
                
                if (rule=="NIDO") {
                      paramname <- paste ("l_", a ,sep="")
                }
                kernel_item_constrain <- paste (kernel_item_constrain, "+", paramname ,sep="")
                
                # add new parameters for the NEW model constrain section
                if (paramname%in%NEW_params_Mplus == FALSE) 
                  NEW_params_Mplus = c(NEW_params_Mplus, paramname)
                
                # add the order restrictions 
                order_item_constrain <- paste(paramname, ">0;", sep="")            
                if (rule!="NIDO") {
                  order_item_constrain <- paste(order_item_constrain, "                         ! main effect for item ", 
                   i, " and attribute ", a, sep="")
                }   
                if (rule=="NIDO") {
                  order_item_constrain <- paste(order_item_constrain, "                         ! main effect for attribute ", a, sep="")
                }  
                
                if (order_item_constrain%in%all_order_constrains == FALSE)                
                    all_order_constrains = c(all_order_constrains, order_item_constrain)
                     
              }   # if     
            }     #for
          }       # if
    
    
          ################### 2-way interactions
          
          if ((anz_att_items > 1) & (max_level_interaction > 1)) {      
              for (j in 1:anzahl_att) {
                  for (k in 2:anzahl_att) {                              
                    if ( (qmatrix[i,j]==1) & (latent_classes[lc,j]==1) & (qmatrix[i,k]==1) & (latent_classes[lc,k]==1) & (j!=k) & (j<k)) {                    
                        # add the 2 way interaction effects
                        paramname <- paste ("l", i, "_2", j , k, sep="")
                        kernel_item_constrain <- paste (kernel_item_constrain, "+", paramname, sep="")
                
                        # add new parameters for the NEW model constrain section
                        if (paramname%in%NEW_params_Mplus == FALSE) 
                           NEW_params_Mplus = c(NEW_params_Mplus, paramname)
                
                        # add the order restrictions 
                        order_item_constrain1 <- paste(paramname, ">-", "l", i, "_1", j, ";", sep="")
                        order_item_constrain1 <- paste(order_item_constrain1, "                   ! 2-way interaction for item ", i, 
                            " and attribute ", j, ",", k, sep="")
                        order_item_constrain2 <- paste(paramname, ">-", "l", i, "_1", k, ";", sep="")
                        order_item_constrain2 <- paste(order_item_constrain2, "                   ! 2-way interaction for item ", i, 
                            " and attribute ", j, ",", k, sep="")      
                
                        if (order_item_constrain1%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain1)                
                        
                        if (order_item_constrain2%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain2)                
                    }   #if                         
                  }     # for
              }         # for
          }             # if
    
    
          ################### 3-way interactions
    
          if ( (anz_att_items > 2) & (max_level_interaction > 2) ) {
              for (j in 1:anzahl_att) {
                for (k in 2:anzahl_att) {
                  for (l in 3:anzahl_att) {                              
                    if ( (qmatrix[i,j]==1) & (latent_classes[lc,j]==1) & (qmatrix[i,k]==1) & (latent_classes[lc,k]==1) & 
                          (qmatrix[i,l]==1) & (latent_classes[lc,l]==1)& (j!=k) & (j!=l) & (k!=l) & (j<k) & (k<l) ) {
                        
                        # add the 3-way interaction effects
                        paramname <- paste ("l", i, "_3", j , k, l, sep="")
                        kernel_item_constrain <- paste (kernel_item_constrain, "+", paramname, sep="")
    
                        # add new parameters for the NEW model constrain section
                        if (paramname%in%NEW_params_Mplus == FALSE) 
                           NEW_params_Mplus = c(NEW_params_Mplus, paramname)
                
                        # add the order restrictions
                        order_item_constrain1 <- paste(paramname, ">-(", "l",i,"_1",j, "+l",i,"_2",j,k, "+l",i,"_2",j,l, ");", sep="")
                        order_item_constrain1 <- paste(order_item_constrain1, "  ! 3-way interaction for item ", i, 
                            " and attributes ", j,",", k,",", l, sep="") 
    
                        order_item_constrain2 <- paste(paramname, ">-(", "l",i,"_1",k, "+l",i,"_2", j,k, "+l", i,"_2",k,l, ");", sep="")
                        order_item_constrain2 <- paste(order_item_constrain2, "  ! 3-way interaction for item ", i, 
                            " and attributes ", j,",", k,",", l, sep="")
                            
                        order_item_constrain3 <- paste(paramname, ">-(", "l",i,"_1",l, "+l",i,"_2",j,l, "+l",i,"_2",k,l, ");", sep="")
                        order_item_constrain3 <- paste(order_item_constrain3, "  ! 3-way interaction for item ", i, 
                            " and attributes ", j,",", k,",", l, sep="")
                
                        if (order_item_constrain1%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain1)                
                        
                        if (order_item_constrain2%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain2) 
                          
                        if (order_item_constrain3%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain3) 
                    
                    }   #if                         
                  }     # for
                }       # for
             }          # for
          }             # if
          
          
          ################### 4-way interactions
          
          if ((anz_att_items > 3) & (max_level_interaction > 3) ) {
              for (j in 1:anzahl_att) {
                for (k in 2:anzahl_att) {
                  for (l in 3:anzahl_att) {                              
                   for (n in 4:anzahl_att) {
                  
                    if ( (qmatrix[i,j]==1) & (latent_classes[lc,j]==1) & (qmatrix[i,k]==1) & (latent_classes[lc,k]==1) & 
                          (qmatrix[i,l]==1) & (latent_classes[lc,l]==1)& (qmatrix[i,n]==1) & (latent_classes[lc,n]==1) &
                          (j!=k) & (j!=l) & (j!=n) & (k!=l) & (k!=n) & (l!=n) & (j<k) & (k<l) & (l<n) ) {
                        
                        # add the 4-way interaction effects
                        paramname <- paste ("l", i, "_4", j , k, l, n, sep="")
                        kernel_item_constrain <- paste (kernel_item_constrain, "+", paramname, sep="")
    
                        # add new parameters for the NEW model constrain section
                        if (paramname%in%NEW_params_Mplus == FALSE) 
                           NEW_params_Mplus = c(NEW_params_Mplus, paramname)
                
                        # add the order restrictions 
                        
                        order_item_constrain1 <- paste(paramname, ">-(", "l",i,"_1",j, "+l",i,"_2",j,k, "+l",i,"_2",j,l,   
                            "+l",i,"_2",j,n, "+l",i,"_3",j,k,l, "+l",i,"_3",j,k,n, "+l",i,"_3",j,l,n, ");", sep="")
                        order_item_constrain1 <- paste(order_item_constrain1, "   ! 4-way interaction for item ", i, 
                            " and attributes ", j, ",", k, "," , l, ",", n, sep="")                    
                        
                        order_item_constrain2 <- paste(paramname, ">-(", "l",i,"_1",k, "+l",i,"_2",j,k, "+l",i,"_2",k,l, 
                            "+l",i,"_2",k,n, "+l",i,"_3",j,k,l, "+l",i,"_3",j,k,n, "+l",i,"_3",k,l,n, ");", sep="")
                        order_item_constrain2 <- paste(order_item_constrain2, "   ! 4-way interaction for item ", i, 
                            " and attributes ", j, ",", k, ",", l, ",", n, sep="")                    
                        
                        order_item_constrain3 <- paste(paramname, ">-(", "l",i,"_1",l, "+l",i,"_2",j,l, "+l",i,"_2",k,l, 
                            "+l",i,"_2",l,n, "+l",i,"_3",j,k,l, "+l",i,"_3",k,l,n, "+l",i,"_3",j,l,n, ");", sep="")
                        order_item_constrain3 <- paste(order_item_constrain3, "   ! 4-way interaction for item ", i, 
                            " and attributes ", j,",", k, ",", l,",", n, sep="")                                            
                            
                        order_item_constrain4 <- paste(paramname, ">-(", "l",i,"_1",n, "+l",i,"_2",j,n, "+l",i,"_2",k,n,   
                            "+l",i,"_2",l,n, "+l",i,"_3",j,k,n, "+l",i,"_3",j,l,n,"+l", i,"_3",k,l,n, ");", sep="")   
                        order_item_constrain4 <- paste(order_item_constrain4, "   ! 4-way interaction for item ", i, 
                            " and attributes ", j, ",", k, ",", l, ",", n, sep="")    
                
                        if (order_item_constrain1%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain1)                
                        
                        if (order_item_constrain2%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain2) 
                          
                        if (order_item_constrain3%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain3) 
                          
                        if (order_item_constrain4%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain4)   
                    
                    }   #if   
                   }    # for                       
                  }     # for
                }       # for
             }          # for
          }             # if
                               
    
          ################### 5-way interactions
    
          if ((anz_att_items > 4) & (max_level_interaction > 4)  ) {
              for (j in 1:anzahl_att) {
                for (k in 2:anzahl_att) {
                  for (l in 3:anzahl_att) {                              
                   for (n in 4:anzahl_att) {
                    for (o in 5:anzahl_att) {
                  
                    if ( ((qmatrix[i,j]==1) & (latent_classes[lc,j]==1) & (qmatrix[i,k]==1) & (latent_classes[lc,k]==1) & 
                          (qmatrix[i,l]==1) & (latent_classes[lc,l]==1)& (qmatrix[i,n]==1) & (latent_classes[lc,n]==1) &
                          (qmatrix[i,o]==1) & (latent_classes[lc,o]==1) &
                          (j!=k) & (j!=l) & (j!=n) & (j!=o) & (k!=l) & (k!=n) & (k!=o) & (l!=n) & (l!=o) &
                          (j<k) & (k<l) & (l<n) & (n<o) ) ) {
                        
                        # add the 5-way interaction effects
                        paramname <- paste ("l", i, "_5", j , k, l, n, o, sep="")
                        kernel_item_constrain <- paste (kernel_item_constrain, "+", paramname, sep="")
    
                        # add new parameters for the NEW model constrain section
                        if (paramname%in%NEW_params_Mplus == FALSE) 
                           NEW_params_Mplus = c(NEW_params_Mplus, paramname)
                
                        # add the order restrictions 
                        
                        order_item_constrain1 <- paste(paramname, ">-(", "l", i, "_1", j, 
                          "+l", i, "_2", j, k, "+l", i, "_2", j, l, "+l", i, "_2", j, n, "+l", i, "_2", j, o,
                          "+l", i,"_3", j,k,l,"+l",i,"_3", j,k,n,"+l", i, "_3", j,k,o,"+l", i,"_3", j,l,n,"+l", i, "_3", j,l,o,"+l", i,"_3", j,n,o,
                          "+l", i,"_4", j,k,l,n, "+l", i,"_4", j,k,l,o, "+l", i,"_4", j,k,n,o, "+l", i,"_4", j,l,n,o, ");", sep="")
                        order_item_constrain1 <- paste(order_item_constrain1, " ! 5-way interaction for item ", i, sep="")
    
                        order_item_constrain2 <- paste(paramname, ">-(", "l", i, "_1", k, 
                          "+l", i, "_2", j, k, "+l", i, "_2", k, l, "+l", i, "_2", k, n, "+l", i, "_2", k, o,
                          "+l", i,"_3", j,k,l,"+l",i,"_3", j,k,n,"+l", i, "_3", j,k,o,"+l", i,"_3", k,l,n,"+l", i, "_3", k,l,o,"+l", i,"_3", k,n,o,
                          "+l", i,"_4", j,k,l,n, "+l", i,"_4", j,k,l,o, "+l", i,"_4", j,k,n,o, "+l", i,"_4", k,l,n,o, ");", sep="")
                        order_item_constrain2 <- paste(order_item_constrain2, " ! 5-way interaction for item ", i, sep="")
                        
                        order_item_constrain3 <- paste(paramname, ">-(", "l", i, "_1", l, 
                          "+l", i, "_2", j,l, "+l", i, "_2", k,l, "+l", i, "_2", l,n, "+l", i, "_2", l,o,
                          "+l", i,"_3", j,k,l,"+l",i,"_3", j,l,n,"+l", i, "_3", j,l,o,"+l", i,"_3", k,l,n,"+l", i, "_3", k,l,o,"+l", i,"_3", l,n,o,
                          "+l", i,"_4", j,k,l,n, "+l", i,"_4", j,k,l,o, "+l", i,"_4", j,l,n,o, "+l", i,"_4", k,l,n,o, ");", sep="")
                        order_item_constrain3 <- paste(order_item_constrain3, " ! 5-way interaction for item ", i, sep="")
           
                        order_item_constrain4 <- paste(paramname, ">-(", "l", i, "_1", n, 
                          "+l", i, "_2", j,n, "+l", i, "_2", k,n, "+l", i, "_2", l,n, "+l", i, "_2", n,o,
                          "+l", i,"_3", j,k,n,"+l",i,"_3", j,l,n,"+l", i,"_3", j,n,o,"+l", i,"_3", k,l,n,"+l", i, "_3", k,n,o, "+l", i,"_3", l,n,o,
                          "+l", i,"_4", j,k,l,n, "+l", i,"_4", j,k,n,o, "+l", i,"_4", j,l,n,o, "+l", i,"_4", k,l,n,o, ");", sep="")
                        order_item_constrain4 <- paste(order_item_constrain4, " ! 5-way interaction for item ", i, sep="")
    
                        order_item_constrain5 <- paste(paramname, ">-(", "l", i, "_1", o, 
                          "+l", i, "_2", j,o, "+l", i, "_2", k,o, "+l", i, "_2", l,o, "+l", i, "_2", n,o,
                          "+l", i,"_3", j,k,o,"+l",i,"_3", j,l,o,"+l", i,"_3", j,n,o,"+l", i,"_3", k,l,o,"+l", i, "_3", k,n,o, "+l", i,"_3", l,n,o,
                          "+l", i,"_4", j,k,l,o, "+l", i,"_4", j,k,n,o, "+l", i,"_4", j,l,n,o, "+l", i,"_4", k,l,n,o, ");", sep="")
                        order_item_constrain5 <- paste(order_item_constrain5, " ! 5-way interaction for item ", i, sep="")                    
                
                        if (order_item_constrain1%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain1)                
                        
                        if (order_item_constrain2%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain2) 
                          
                        if (order_item_constrain3%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain3) 
                          
                        if (order_item_constrain4%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain4) 
                          
                        if (order_item_constrain5%in%all_order_constrains == FALSE) 
                          all_order_constrains = c(all_order_constrains, order_item_constrain5)    
                    
                    }   #if   
                    }   # for
                   }    # for                       
                  }     # for
                }       # for
             }          # for
          }             # if
                                        
          
          # If the kernel term is already known, find the right position in the array with known constrains
          # and get the name of the kernel
          if (kernel_item_constrain%in%kernel_item_without_names == TRUE)  {                                
               index <- which(kernel_item_without_names == kernel_item_constrain)          
               # get the name of the saved kernal
               substrings <- unlist(strsplit(as.character(all_kernel_constrains[index]), "=", fixed = TRUE))
               kernel_item_name <- substrings[1]
               
               kernel_item_constrain <- paste(kernel_item_name, "=", kernel_item_constrain, sep="")                                            
          }  
          
          # If the kernal term is new, create a new name and save the new kernal 
          if (kernel_item_constrain%in%all_kernel_constrains == FALSE)  {                             
               m <- m +1         # iterator for the names of the kernals for the current item defined above
               kernel_item_name <- paste("t", i, "_", m, sep="")
               # save kernal without name
               kernel_item_without_names <- c(kernel_item_without_names, kernel_item_constrain)
               kernel_item_constrain <- paste(kernel_item_name, "=", kernel_item_constrain, sep="")
               # save kernal with the name
               all_kernel_constrains = c(all_kernel_constrains, kernel_item_constrain)                              
          } 
           
          # needed only for test purposes
          kernel_namen[i,lc] <- kernel_item_name           
          kernel_item_constrains[i,lc] <- kernel_item_constrain                       
    
      }      #### loop over latent classes
    }       #### loop over items
    
    # delete the first empty string from each of the vectors
    all_kernel_constrains <- all_kernel_constrains[2:length(all_kernel_constrains)]
    all_order_constrains <- all_order_constrains[2:length(all_order_constrains)]
    NEW_params_Mplus <-  NEW_params_Mplus[2:length(NEW_params_Mplus)]
    
    # Multiply all kernels with -1 because of differences in Mplus program
    all_kernel_constrains <- gsub("=", "=-(", all_kernel_constrains)
    all_kernel_constrains <- paste(all_kernel_constrains, ");", sep="")
    
    
    #--------------------------------------------------------------------------------------------------
    #------- Prepare data for the Mplus-Syntax --------------------------------------------------------
    
    #### Produce a matrix with all variable names incl. PreAddVars and PostAddVars
    # !!! Problem: lines can get longer than 80 characters
    # Solution: write variable names in multiple lines
    varnames_all = PreAddVars
    varnames_multiple_lines_all <- matrix(data=NA, 52, 1)
    num_chars_line <- 1
    num_chars_count <- nchar(varnames_all)
    for (r in 1:anzahl_items) {            
            varnames_all <- paste(varnames_all, " ", rownames(qmatrix)[r], sep="")
            num_chars_count <- num_chars_count + nchar(rownames(qmatrix)[r])+1
            
            if ((num_chars_count > 55) & (num_chars_line<length(varnames_multiple_lines_all[,1])) ) {
                num_chars_count <-0               
                varnames_multiple_lines_all[num_chars_line,1]<- varnames_all
                num_chars_line <- num_chars_line + 1
                varnames_all=""
            }
    }    
    
    if ((varnames_all=="") & (PostAddVars=="") ) 
        { varnames_multiple_lines_all[num_chars_line-1,1]<- paste(varnames_multiple_lines_all[num_chars_line-1,1], ";", sep="") }     
        
    if ((varnames_all!="") & (PostAddVars=="") ) 
        { varnames_multiple_lines_all[num_chars_line,1]<- paste(varnames_all, ";", sep="") }
        
    if ((varnames_all!="") & (PostAddVars!="") ) 
        { varnames_multiple_lines_all[num_chars_line,1]<- varnames_all }        
        
    if (PostAddVars!="") { varnames_multiple_lines_all[num_chars_line+1,1]<- paste(PostAddVars, ";", sep="") }

    #### Produce a matrix with only relevant variable names
    # !!! Problem: lines can get longer than 80 characters
    # Solution: write variable names in multiple lines
    varnames = ""
    varnames_multiple_lines <- matrix(data=NA, 50, 1)
    num_chars_line <- 1
    num_chars_count <-0
    for (r in 1:anzahl_items) {            
            varnames <- paste(varnames, " ", rownames(qmatrix)[r], sep="")
            num_chars_count <- num_chars_count + nchar(rownames(qmatrix)[r])+1
            
            if ((num_chars_count > 55) & (num_chars_line<length(varnames_multiple_lines[,1])) ) {
                num_chars_count <-0               
                varnames_multiple_lines[num_chars_line,1]<- varnames
                num_chars_line <- num_chars_line + 1
                varnames=""
            }
    }    
    if (varnames!="") { varnames_multiple_lines[num_chars_line,1]<- paste(varnames, ";", sep="") }  
    if (varnames=="") { varnames_multiple_lines[num_chars_line-1,1]<- paste(varnames_multiple_lines[num_chars_line-1,1], ";", sep="") }             
    
    #### Write the complete pattern for each latent class in a new matrix
    pattern_latent_classes = matrix(data=NA, max_att, 2)
    for (lc in 1: max_att) {
        pattern_latent_classes[lc,1] <- lc
        pattern_latent_classes[lc,2] <- ""
        for (a in 1:anzahl_att) {
            pattern_latent_classes[lc,2] <- paste(pattern_latent_classes[lc,2], latent_classes[lc,a], sep="")
        }
    }
    colnames(pattern_latent_classes) = c("Latent class", "attribute pattern")
    
    #### Write the complete pattern for each item in a new matrix
    pattern_item <- matrix(data=NA, anzahl_items, 2)
    for (i in 1:anzahl_items) {
        pattern_item[i,1] <- i
        pattern_item[i,2] <- ""
        for (a in 1:anzahl_att) {
            pattern_item[i,2] <- paste(pattern_item[i,2], qmatrix[i,a], sep="")
        }
    }
    colnames(pattern_item) = c("Item", "Qmatrix pattern")
    
    #### Find the item number for all item kernals 
    item_kernal <- matrix(data=NA, length(all_kernel_constrains), 2)
    for (m in 1:(length(all_kernel_constrains))) {
      
        # get the number of the item
        substrings <- unlist(strsplit(as.character(all_kernel_constrains[m]), "_", fixed = TRUE))
        substrings_neu <- unlist(strsplit(as.character(substrings[1]), "t", fixed = TRUE))         
           
        item_kernal[m,1] <- substrings_neu[2]         
        item_kernal[m,2] <- all_kernel_constrains[m]
    } 
    
    if (rule!="NIDO") {
        #### Find the item number for all order constains 
        item_order <- matrix(data=NA, length(all_order_constrains), 2)
        for (m in 1:(length(all_order_constrains))) {
             
             # get the number of the item
             substrings <- unlist(strsplit(as.character(all_order_constrains[m]), "_", fixed = TRUE))
             substrings_neu <- unlist(strsplit(as.character(substrings[1]), "l", fixed = TRUE))         
             
             item_order[m,1] <- substrings_neu[2]         
             item_order[m,2] <- all_order_constrains[m]
        } 
    
        #### Find the item number for all parameters in NEW-vector
        item_new_param <- matrix(data=NA, length(NEW_params_Mplus), 2)
        for (n in 1:(length(NEW_params_Mplus))) {
            
             # get the number of the item
             substrings <- unlist(strsplit(as.character(NEW_params_Mplus[n]), "_", fixed = TRUE))
             substrings_neu <- unlist(strsplit(as.character(substrings[1]), "l", fixed = TRUE)) 
             
             item_new_param[n,1] <- substrings_neu[2]
             item_new_param[n,2] <- NEW_params_Mplus[n]
        } 
    }
    
    #------------------------------------------------------------------------------------------
    #------- Creating Syntax for Mplus --------------------------------------------------------
    
    sink(input_file)
    
    ### TITLE section
    print("TITLE:", quote =FALSE)
    print(paste("! DCM for ", filename_ds, " with ", anzahl_att, " attributes and ", anzahl_items, " items", sep=""), quote =FALSE)
    if (rule == "FULL") {
        print(paste("! Saturated structural model with ", max_level_interaction, " levels of interactions", sep=""), quote =FALSE)
    }
    if (rule != "FULL") {
        print(paste("! LCDM model with ", rule, "-rule", sep=""), quote =FALSE)
    }       
    print ("",quote =FALSE)    
    
    ### DATA section
    print("DATA:  ! location of free format data file (in syntax file folder);", quote =FALSE)
    print(paste("FILE IS ", filename_ds, ";", sep=""), quote =FALSE)
    print ("",quote =FALSE)
    
    ### VARIABLE section
    print("VARIABLE:",quote =FALSE)
    
    # NAMES ARE    
    print( "NAMES ARE " ,quote =FALSE)  
    # print multiple lines if the matrix 'varnames_multiple_lines_all' is not empty     
    for (i in 1:length(varnames_multiple_lines_all[,1]) ) {    
        if (is.na(varnames_multiple_lines_all[i,1])==F) {
            print(varnames_multiple_lines_all[i,1],quote =FALSE)      
        }
    }
    # print multiple lines if the matrix 'varnames_multiple_lines' is not empty
    print ("",quote =FALSE)
        
    # USEVARIABLE ARE
    print( "USEVARIABLE ARE" ,quote =FALSE)
    # print multiple lines if the matrix 'varnames_multiple_lines' is not empty    
    for (i in 1:length(varnames_multiple_lines[,1]) ) {    
        if (is.na(varnames_multiple_lines[i,1])==F) {
            print(varnames_multiple_lines[i,1],quote =FALSE)      
        }
    }
    print ("",quote =FALSE)    
    
    # CATEGORICAL ARE
    print( "CATEGORICAL ARE" ,quote =FALSE)
    # print multiple lines if the matrix 'varnames_multiple_lines' is not empty    
    for (i in 1:length(varnames_multiple_lines[,1]) ) {    
        if (is.na(varnames_multiple_lines[i,1])==F) {
            print(varnames_multiple_lines[i,1],quote =FALSE)      
        }
    }
    print ("",quote =FALSE) 
    
    # if requested MISSING ARE
    if (miss==T) {
          print( paste("MISSING ARE ALL (", miss_char, ");" , sep=""), quote =FALSE)
          print ("",quote =FALSE)
    } 
    
    # CLASSES
    print( paste("CLASSES = c(", max_att, ");    ! ",max_att, " possible attribute patterns for model with ", anzahl_att, " attributes", sep="") ,quote =FALSE)
    print ("",quote =FALSE)
    
    ### ANALYSIS section
    print("ANALYSIS:", quote =FALSE)
    print("TYPE=MIXTURE;     ! estimates latent classes;",quote =FALSE)
    print("STARTS=0;         ! turn off multiple random start feature (disabled anyway);", quote =FALSE)        
    print ("",quote =FALSE)
    
    ### MODEL section
    print("MODEL:", quote =FALSE)
    print("%OVERALL%", quote =FALSE)
    print("", quote =FALSE)
    
    for (lc in 1:max_att) {
        print(paste("%c#", lc, "%     ! for latent class ", lc, "=(", pattern_latent_classes[lc,2], ")", sep="") , quote =FALSE) 
        
        for (i in 1:anzahl_items) {
          print(paste("[", rownames(qmatrix)[i], "$1]  (", kernel_namen[i,lc], ");     ! threshold for item ", i, sep=""), quote =FALSE)            
        }       
    }    
    print("", quote =FALSE)
    
    ### MODEL CONSTRAINT section
    print("MODEL CONSTRAINT:          ! used to define LCDM parameters and constrains", quote =FALSE)
    print("! NOTE: Mplus uses P(X=0) rather than P(X=1) so terms must be multiplied by -1", quote =FALSE)
    print("", quote =FALSE)
    
    # If the rule is NIDO print all attribute specific order constrains right after the 
    # NEW statement
    if (rule=="NIDO") {
        new_terms = ""
        for (r in 1:length(NEW_params_Mplus)) {                          
            new_terms <- paste(new_terms, " ", NEW_params_Mplus[r], sep="")                                      
        }
        print(paste("NEW (", new_terms ,");", sep=""), quote =FALSE)        
                       
        for (r in 1:length(all_order_constrains)) {                          
            print(all_order_constrains[r], quote=FALSE)                    
        }        
    }     
    
    #### Print NEW expression, kernal expressions and order restriction for each item
    for (i in (1:anzahl_items)) {
    
        print(paste("! ITEM ", i, ":", sep=""), quote =FALSE)
        print(paste ("! Q-matrix Entry ", pattern_item[i,2], sep=""), quote =FALSE)                           
         
       if (rule!="NIDO") {        
            #### Print the NEW command by controlling the length of the line      
            # !!! Problem: lines can get longer than 90 characters
            # Solution: write new terms in multiple lines
            new_terms = "NEW ("
            new_terms_multiple_lines <- matrix(data=NA, 8, 1)
            num_chars_line <- 1
            num_chars_count <-nchar(new_terms)
            
            for (r in 1:length(item_new_param[,1])) {                    
              if (item_new_param[r,1]==i) {            
                    
                    # write the new terms in a new matrix      
                    new_terms <- paste(new_terms, " ", item_new_param[r,2], sep="")
                    num_chars_count <- num_chars_count + nchar(item_new_param[r,2])+1
                            
                    if ((num_chars_count > 70) & (num_chars_line < length(varnames_multiple_lines[,1])) ) {
                            num_chars_count <-0               
                            new_terms_multiple_lines[num_chars_line,1]<- new_terms
                            num_chars_line <- num_chars_line + 1
                            new_terms=""
                    }                           
              }  
            } # for    
            if (new_terms!="") { new_terms_multiple_lines[num_chars_line,1]<- paste(new_terms, ");", sep="") }                
                    
            # print the matrix "new_terms_multiple_lines"
            for (w in 1:length(new_terms_multiple_lines[,1]) ) {    
                if (is.na(new_terms_multiple_lines[w,1])==F) {
                      print(new_terms_multiple_lines[w,1],quote =FALSE)      
                }
            }        
        } # rule!=NIDO
        
        #### print the kernal constrain for the current item
        # !!! Problem: lines can get longer than 90 characters
        # Solution: write new terms in multiple lines
        for (r in 1:length(item_kernal[,1])) {    
            if ( (item_kernal[r,1]==i) & (nchar(item_kernal[r,2])<89) ) { 
               print(item_kernal[r,2],quote=FALSE)            
            }
                        
            if ( (item_kernal[r,1]==i) & (nchar(item_kernal[r,2])>89) ) {                         
                                          
               # write the new terms in a new matrix
               substrings <- unlist(strsplit(as.character(item_kernal[r,2]), "+", fixed = TRUE))
               new_kernal = substrings[1]
               new_kernal_multiple_lines <- matrix(data=NA, 8, 1)
               num_chars_line <- 1
               num_chars_count <-nchar(new_kernal)                              
               
               for (z in 2:length(substrings)) {                        
                  new_kernal <- paste(new_kernal, "+", substrings[z], sep="")
                  num_chars_count <- num_chars_count + nchar(substrings[z])+1
                           
                  if ((num_chars_count > 70) & (num_chars_line < length(new_kernal_multiple_lines[,1])) ) {                                       
                        new_kernal_multiple_lines[num_chars_line,1]<- new_kernal
                        num_chars_count <-0
                        num_chars_line <- num_chars_line + 1
                        new_kernal=""
                  }  
               }    
               if (new_kernal!="") { new_kernal_multiple_lines[num_chars_line,1]<- new_kernal }                           
               
               # print the shortened lines
               for (w in 1:length(new_kernal_multiple_lines[,1]) ) {    
                  if (is.na(new_kernal_multiple_lines[w,1])==F) {
                      print(new_kernal_multiple_lines[w,1],quote =FALSE)      
                  }
               }                                                         
            }                                    
        } 
        
        if (rule!="NIDO" & orderConstrains==TRUE) {
            # print the order restrictions for the current item
            print("! Order constrains", quote =FALSE)
            for (r in 1:length(item_order[,1])) {    
                if (item_order[r,1]==i) { 
                  print(item_order[r,2],quote=FALSE)            
                }
            }                    
            print("", quote =FALSE)    
        }
    }     # for i
    
    ### OUTPUT section
    print("", quote =FALSE)
    print("OUTPUT:", quote =FALSE)
    print("   TECH10;", quote =FALSE)
    
    newname <- gsub(".inp", "", input_file)
    ### SAVEDATA section
    print("", quote =FALSE)
    print("SAVEDATA:", quote =FALSE)
    print("   FORMAT is f10.5;                 ! format of output file", quote =FALSE)
    print(paste("   FILE is ", newname, ".cprob;         ! specify file name and location", sep=""), quote =FALSE)
    print("   SAVE = CPROBABILITIES;           ! specify the information to save", quote =FALSE) 
    # print(paste("   ESTIMATES ARE ", newname, ".est;     ! save the estimates", sep="" ), quote =FALSE) 
    print(paste("   RESULTS ARE ", newname, ".res;       ! save the results", sep="" ), quote =FALSE) 
        
    sink()
    
    ### Cut of the line numbers printed to the output
    mplus_input <- readLines(input_file)
    # meine Variante
    mplus_input <- gsub("[//[]1]", "", mplus_input) 
    # Wichts variante
    # mplus_input <- gsub("^\\s*\\[1\\]\\s*", "", mplus_input) 
    writeLines(mplus_input, input_file)

}
############## End of function generate.mplus.input ####################################



