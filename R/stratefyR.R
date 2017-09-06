#######################################
###       StratefyR function       ####
#######################################
##############################################
### Calculate CSR ecological strategies using
### Pierce et al. 2017 methods
### (Functional Ecology 31(2): 444-457) ###
##############################################
##' CSR ecology strategies calculator
##' 
##' This function...
##'
##' @param sp species names or morphospecies codes
##' @param la leaf area (mm^2) 
##' @param lfw leaf fresh weight (mg)
##' @param ldw leaf dry weight (mg)
##' @param data a data frame with species name, LA, LFW, LDW  values for each species. Names should be 'sp', 'la', 'lfw' e 'ldw'.
##' @return returns a data frame with C, S and R index strategies for each species.
##' @author Renan Parmigiani \and Alexandre Adalardo de Oliveira \email{ecovirtualpackage@@gmail.com} 
##' @keywords ecology strategies, competition, stress-tolerant , ruderal, plants 
## @import tcltk
## @importFrom grDevices colorRamp dev.new rainbow rgb
## @importFrom graphics abline axis curve grid image layout legend lines matplot mtext par plot points polygon segments text title
##' @examples
##' 
##' sp <- c("Molopospermum p.", "Kalmia p.", "Arabidopsis t.", "Pteridium a.")    
##' la <-  c(369615.7, 11.8,55.7,36061.2)                                         
##' lfw <- c(84009.6,4.5,12.3,6988.0)                                             
##' ldw <- c(21187.8,1.8,1.6,2480.0)                                               
##' traitsp <- data.frame(sp, la, lfw, ldw)
##' stratefyR(data = traitsp)     
##' 
##' 
##' @export stratefyR         
stratefyR <- function(sp = NULL, la = NULL, lfw = NULL, ldw = NULL, data = NULL)
{
    if(!is.null(data))
    {
        nsp = nrow(data)
        names(data) <- tolower(names(data))
        la = data$la
        lfw = data$lfw
        ldw = data$ldw
        if(sum(names(data)=="sp")< 1)
        {
            sp = paste("sp", 1:nsp, sep="_")
        }
        else
        {
            sp = data$sp
        }
            
    }
    
if((is.null(la) | is.null(lfw) | is.null(ldw)) )
{
        stop("at least one of the tree main traits (LA, LFW, LDW) are missing")
        
}
    ## Calculated traits
    sla <- la/ldw                      # calculating SLA [L]
    suc.ind <- (lfw - ldw)/(la/10)     # leaf succulence index [AC]
    vfsuc <- suc.ind > 5
    ldmc <- (ldw*100)/lfw
    ldmc[vfsuc] <- (100-((ldw[vfsuc]*100)/lfw[vfsuc]))
    la.sqrt <- sqrt(la/894205)*100       #função de ligação do valor de LA, esse valor colocado na função (894205), faz referência ao valor máximo encontrado no banco de dados usados, esse valor pode ser usado considerando que o trabalho dele contempla toda amplitude das áreas foliares das plantas do mundo (n= 3068) [AE]
    ldmc.logit <- log((ldmc/100)/(1-(ldmc/100)))    #função de ligação LDMC [AF]
    sla.log <- log(sla)                            #função de ligação SLA [AG]
    pca2 <- -0.8678 + 1.6464 * la.sqrt              #modelo que ajusta o LA ao PCA2 [AH]
    #pca1.pos <- 1.3369 + 1.0019e-05 * (1-exp(-2.2303e-12 * ldmc.logit)) + 4.5835*(1-exp(-0.2328 * ldmc.logit))   #modelo que ajusta o LDMC à parte positiva do PCA1 [AI]
    pca1.pos <- 1.3369+0.000010019*(1-exp(-0.0000000000022303*ldmc.logit)) + 4.5835*(1 - exp(-0.2328 * ldmc.logit)) 
    pca1.neg <- -57.5924 + 62.6802 * exp(-0.0288 * sla.log)    #modelo que ajusta o SLA à parte negativa do PCA1 [AJ]
##########################################
### Ajusting for minimal or maximum values
##########################################
    min.C <- 0                        #menor valor do eixo C segundo os dados do stratefy (Pierce et al., 2017) [AK]
    min.S <- -0.756451214853076       #menor valor do eixo S segundo os dados do stratefy (Pierce et al., 2017) [AL]
    min.R <- -11.3467682227961        #menor valor do eixo R segundo os dados do stratefy (Pierce et al., 2017) [AM]
    neg.outl.C <- pca2 ## pode usar direto o objeto pca2
    neg.outl.C[pca2 < min.C] <-  min.C              #caso o valor calculado (pca2) seja inferior ao valor mínimo dos dados do stratefy (Pierce et al., 2017) usar o valor mínimo [AN]
    neg.outl.S <- pca1.pos
    neg.outl.S[pca1.pos < min.S]  <- min.S                  #ajuste do valor calculado do eixo S para fi
    neg.outl.R <- pca1.neg
    neg.outl.R[pca1.neg < min.R]  <- min.R                  #ajuste do valor calculado do eixo S para fi
################
    max.C <- 57.3756711966087           #maior valor do eixo C  segundo os dados do stratefy (Pierce et al., 2017) [AQ]
    max.S <- 5.79158377609218           #IDEM para o eixo S [AR] 
    max.R <- 1.10795515716546
    neg.outl.C[pca2 > max.C] <-  max.C              #caso o valor calculado (pca2) seja inferior ao valor mínimo dos dados do stratefy (Pierce et al., 2017) usar o valor mínimo [AN]
    neg.outl.S[pca1.pos > max.S]  <- max.S                  #ajuste do valor calculado do eixo S para fi
    neg.outl.R[pca1.neg > max.R]  <- max.R   
    min.pos.C <- abs(min.C)             #tornando os valores de C mínimos positivos para a amplitude se dar dentro de um intervalo positivo [AW]
    min.pos.S <- abs(min.S)             #IDEM para o eico S [AX]
    min.pos.R <- abs(min.R)             #IDEM para o eico R [AY]
    valor.C <- min.pos.C + neg.outl.C   #Cálculo do valor de C corrigido somando a parte negativa [AZ]
    valor.S <- min.pos.S + neg.outl.S   #Cálculo do valor de S corrigido somando a parte negativa [BA]
    valor.R <- min.pos.R + neg.outl.R   #Cálculo do valor de R corrigido somando a parte negativa [BB]
    range.C <- max.C + abs(min.pos.C)   #Intervalo de valores entre o mínimo e o máximo do eixo C [BC]
    range.S <- max.S + abs(min.pos.S)   #Intervalo de valores entre o mínimo e o máximo do eixo S [BD]
    range.R <- max.R + abs(min.pos.R)   #Intervalo de valores entre o mínimo e o máximo do eixo R [BE]
    prop.C  <- (valor.C/range.C)*100        #Proporção da varição do eixo C [BF]
    prop.S  <- (valor.S/range.S)*100        #Proporção da varição do eixo S [BG]
    prop.R  <- 100-((valor.R/range.R)*100)  #Proporção da varição do eixo R [BH]
    conv.porc <- 100/(apply(cbind(prop.C, prop.S, prop.R), 1, sum))    #valor que converte a proporção em porcentagem [BI]
    C <- prop.C * conv.porc   #Valor do eixo C em porcentagem [M]
    S <- prop.S * conv.porc   #Valor do eixo S em porcentagem [N]
    R <- prop.R * conv.porc   #Valor do eixo R em porcentagem [O]
  return(data.frame(C,S,R))
}
######################END############################
