possibilities <- function(n,states=3){
        res <- NULL
        for (i in n:0){
                a <- i
                if (i==n){
                        b <- 0
                        c <- 0
                        res <- rbind(res,c(a,b,c))
                }
                else{
                        b <- n-a
                        while(b>=0){
                                c <- (n-a-b)
                                res <- rbind(res,c(a,b,c))
                                b <- b-1
                        }
                }
           
        }
        res
}

trinom <- function(n, state, prob){
        
        factorial(n)/(factorial(state[,1])*factorial(state[,2])*factorial(state[,3]))*
                (prob[1]^state[,1])*(prob[2]^state[,2])*(prob[3]^state[,3])
        
}

Outage_Calc <- function(Total_MW,Units,prob_out,prob_par,PDR,MTTR_O,MTTR_P){
        Av_MW <- Total_MW /  Units 
        P_MW <- Av_MW * (1-PDR)
        
        prob_avail <- (1-prob_out)*(1-prob_par)
        prob_pout <- (1-prob_out)*prob_par
        
        unit_config <- possibilities(Units)
        unit_rating <- c(Av_MW,P_MW,0)
        prob_rating <- c(prob_avail,prob_pout,prob_out)
        
        MTTR <- unit_config[,2:3]%*%c(MTTR_O,MTTR_P)/(unit_config[,2]+unit_config[,3])
        MTTR[is.na(MTTR)] <- 0
        OutFac <-  round((Total_MW - unit_config%*%unit_rating)/Total_MW,3)
        FOR <- round(trinom(Units,unit_config,prob_rating),6)
        Final <- cbind(unit_config,FOR,OutFac,MTTR )
        colnames(Final) <- c("Avail","Partial Out","Full Out", "FOR","Outage Factor","MTTR")
        Final
        
        }



Murray <- Outage_Calc(1500,14,0.0082217,0.0000333,0.119697,33.45*2,10.75*2)
Hallet <- Outage_Calc(220,12,0.0066015,0.0008687,0.3189,20.45*2,5.9*2)

Murray[Murray[,4]>0,]
