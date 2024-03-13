#' Variations
#'
#' @description
#' The Variations function provides a visualization about the effect of different parameters on the Rank-abundance curve.
#'
#'
#' @usage Variaciones(sp, r, mc, mp, z, y, x)
#'
#'
#'
#' @param sp Species number of the simulator
#' @param r Number of resources in the simulator
#' @param mc Media of consumption of the species
#' @param mp Media of production of the species
#' @param z The parameter that will change in next iteration, 1= sp, 2= r, 3= mc, 4= mp.
#' @param y Value that will be add in each iteration
#' @param x Number of iterations
#'
#' @author Erik Martinez
#'
#' @example C:/Users/DELL/Desktop/GitHub/Funcion_Rank_Abundance_Variaciones/FuncionR/R/Example_function.R
#'
#'
#'
#'
#' @returns A data frame that contains the A,a, b values of the different Rank-Abundance curves
#'
#'



Variaciones <- function(sp,r,mc,mp,z,y,x){

  library(miaSim)
  library(DGBD)
  library(ggplot2)
  library(RADanalysis)


  n_species <- sp
  n_resources <-r
  library(miaSim)
  library(DGBD)
  library(ggplot2)
  library(RADanalysis)
  valoresfinales_Aab <- c()
  parametros_simulador <- data.frame()
  Base_conjunta_de_variaciones <- data.frame()


  for (i in 1:x) {



    CRsims <- list()
    for(i in 1:10){ #aquii el numero de iteraciones


      matE <- randomE(
        n_species = n_species, n_resources = n_resources,
        mean_consumption = mc, mean_production = mp, maintenance = 0.4
      )

      resources_c <- rep(100, n_resources)
      CRsims[[i]] <- simulateConsumerResource(
        n_species = n_species,
        n_resources = n_resources, names_species = letters[seq_len(n_species)],
        names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = matE,
        x0 = rep(0.001, n_species), resources = resources_c,
        growth_rates = runif(n_species),
        error_variance = 0.001,
        t_end = 15000,
        t_step = 1,
        # t_external_events = 5000, #evento externo en que momento
        #t_external_durations = 300, #Cuanto va a durar ese evento
        sigma_external = 0.05,
        stochastic = T,
        norm = T
      )
    }

    ab_tables_div<-function(ab_tables_list, diversity_type){
      if(diversity_type == "shannon"){
        div_list<-list()
        for(j in 1:length(ab_tables_list)){
          div_table<-c()
          for(i in 1:ncol(ab_tables_list[[j]])){
            div_table[i]<-Div_Shannon(ab_tables_list[[j]][, i])
          }
          div_list[[j]]<-div_table
        }
      } else
        if(diversity_type == "simpson"){
          div_list<-list()
          for(j in 1:length(ab_tables_list)){
            div_table<-c()
            for(i in 1:ncol(ab_tables_list[[j]])){
              div_table[i]<-Dom_Simpson(ab_tables_list[[j]][, i])
            }
            div_list[[j]]<-div_table
          }
        } else
          if(diversity_type == "pielou"){
            div_list<-list()
            for(j in 1:length(ab_tables_list)){
              div_table<-c()
              for(i in 1:ncol(ab_tables_list[[j]])){
                div_table[i]<-Eq_Pielou(ab_tables_list[[j]][, i])
              }
              div_list[[j]]<-div_table
            }
          } else
            if(diversity_type == "ginisimpson"){
              div_list<-list()
              for(j in 1:length(ab_tables_list)){
                div_table<-c()
                for(i in 1:ncol(ab_tables_list[[j]])){
                  div_table[i]<-1-Dom_Simpson(ab_tables_list[[j]][, i])
                }
                div_list[[j]]<-div_table
              }
            }
      return(div_list)
    }


    # Shannon diversity
    Div_Shannon <- function(abundancias_ab){
      abs_rel <- abundancias_ab/sum(abundancias_ab)
      Shannon <- -sum(abs_rel*log(abs_rel))
      return(Shannon)
    }

    # Simpson dominance
    Dom_Simpson <- function(abundancias_ab){
      abs_rel <- abundancias_ab/sum(abundancias_ab)
      Simpson <- sum(abs_rel^2)
      return(Simpson)
    }

    # Pielou evenness
    Eq_Pielou <- function(abundancias_ab){
      abs_rel <- abundancias_ab/sum(abundancias_ab)
      Shannon <- -sum(abs_rel*log(abs_rel))
      Pielou <- Shannon/log(length(abundancias_ab))
      return(Pielou)
    }

    # Extract abundance tables
    CRsims_T<-list()
    for(i in 1:length(CRsims)){
      CRsims_T[[i]]<-as.matrix(assay(CRsims[[i]]))
    }

    # Separate the data
    sep_tables<-function(ab_tables, start_p, end_p){
      sep_tables<-list()
      for(i in 1:length(ab_tables)){
        sep_tables[[i]]<-ab_tables[[i]][,start_p:end_p]
      }
      return(sep_tables)
    }
    ########################################
    basal_tables <- sep_tables(CRsims_T, 300, 1000)

    ###########################################
    # Pielou diversity
    pielou_basal<-unlist(ab_tables_div(basal_tables, "pielou"))


    # Simpson dominance
    simpson_basal<-unlist(ab_tables_div(basal_tables, "simpson"))


    df_diversity<-data.frame(
      Period = rep(1, length(pielou_basal)),
      Pielou =  pielou_basal,
      Simpson = simpson_basal
    )

    # Rank abundance analysis
    rank_abs<-function(ab_tables, start_p, end_p){
      rank_abs<-list()
      for(j in 1:length(ab_tables)){
        sum_abs<-c()
        for(i in 1:nrow(ab_tables[[1]])){
          sum_abs[i]<-sum(ab_tables[[j]][,start_p:end_p][i,])
        }
        rank_abs[[j]]<-sort(sum_abs/nrow(ab_tables[[1]]), decreasing = T)
      }
      rank_abs_mat<-matrix(unlist(rank_abs), nrow=length(ab_tables),
                           ncol=nrow(ab_tables[[1]]), byrow=TRUE)
      return(rank_abs_mat)
    }

    basal_rank <- rank_abs(CRsims_T, 300, 1000)

    #print(basal_rank)


    basalrank_prueba <- c()
    ini <- 1

    for (i in 1:n_species) {

      Prom <- sum(basal_rank[,ini]/10 )

      basalrank_prueba <- rbind(basalrank_prueba,Prom)

      ini <- ini+1
      print(Prom)
    }



    library(DGBD)


    valores_Aab <- BC_multiple(basalrank_prueba)

    parametros_Aab <- valores_Aab[[1]][[2]]

    parametros_simulador <- cbind(parametros_Aab[1], parametros_Aab[2], parametros_Aab[3], sp, r, mc, mp)

    Base_conjunta_de_variaciones <- rbind(Base_conjunta_de_variaciones, parametros_simulador)




    if(z == 1){
      sp<- sp+y
    }else if (z==2){
      r<- r+y
    }else if (z==3){
      mc<- mc+y
    }else if(z==4){
      mp<- mp+y
    }else{
      print("lee las instrucciones")
    }


    print(i)
  } #corchete del for inicial

  colnames(Base_conjunta_de_variaciones) <- c("A","a", "b", "sp","r","mc","mp")
  print(Base_conjunta_de_variaciones)

}




