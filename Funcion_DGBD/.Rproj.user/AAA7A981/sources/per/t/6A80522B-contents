#sp <- numero de especies
#r <- numero de recursos
#mc <- media de consumo MC
#mp <- media de produccion MP
#z <- valor a variar 1=sp, 2=R, 3=MC, 4=MP
#y <- la variacion que se quiere agregar
#x <- numero de veces que quieres repetirlo
devtools::install_github("Fa-Moe/DGBD", force = T)


Variaciones <- function(sp,r,mc,mp,z,y,x){

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

  basal_tables <- sep_tables(CRsims_T, 300, 1000)
  
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
  
basalrank_prueba <- c()

ini <- 1

for (i in 1:n_species) {
  
  Prom <- sum(basal_rank[,ini])/10
  
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



} #corchete del for inicial

  colnames(Base_conjunta_de_variaciones) <- c("A","a", "b", "sp","r","mc","mp")
print(Base_conjunta_de_variaciones)



}


#sp <- numero de especies
#r <- numero de recursos
#mc <- media de consumo MC
#mp <- media de produccion MP
#z <- valor a variar 1=sp, 2=R, 3=MC, 4=MP
#y <- la variacion que se quiere agregar
#x <- numero de veces que quieres repetirlo

#Variaciones(sp,r_,mc,mp,z,y,x) 


#ESPECIES 25, 30, 40 LAP, 50, 75, 100, 200, 300
#RECURSOS 50, 75, 120, 250
#MC   25, 75, 125, 175, 225
#MP 50, 100, 150, 200, 250

data_example <- Variaciones(20,20,25,25,4,2,3)
saveRDS(data_example,"03_result/example4.rds")



#Solo 
#GENERAR UN REPORTE DE DOCUMENTADO DE LA FUNCION 
#WARNINGS
#REPORTE DE LA FUNCION



ejemplo1 <- readRDS("03_result/example1.rds")
ejemplo1
