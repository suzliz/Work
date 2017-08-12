## Use to expand possible operating states for n units of equal size
trinomial_exp <- function(n){
        res <- NULL
        ## a = units available
        ## b = units partially available
        ## c = units on full outage
        for (i in n:0){ 
                ## work down from n units available to zero
                a <- i
                if (i==n){
                        b <- 0
                        c <- 0
                        #Add configuration to results matrix
                        res <- rbind(res,c(a,b,c))
                }
                else{
                        ##Calculate maximum partially available units        
                        b <- n-a
                        while(b>=0){
                                ## Calculate units on full outage
                                ## Ensures sum of units in each state = total units
                                c <- (n-a-b)
                                ## Add configuration to results matrix
                                res <- rbind(res,c(a,b,c))
                                ## Lower units on partial outage by 1
                                ## This increases full outage units by 1
                                b <- b-1
                        }
                }
                
        }
        colnames(res) <- c("Available","Partial","Full_Out")
        res
}

## Use to expand possible operating states for up to 4 units
full_exp <- function(num){
        ## Generate matrix of states -all values set to zero initially
        ## Sorry - Horrible set of nested loops
        
        ## Initialise results matrix
        ## There are 3 ^ n combinations of operating states (rows)
        ## There are 3 * n columns (3 possibilities for each of n units)
        res <- matrix(nrow=3^num,ncol=3*num,0)
        
        ##row counter
        n <- 1 
        ## Calculate the three possible outage states for each unit
        ## Available, Partial Outage and Full Outage
        for(i in 1:3){
                if (num==1){ ## One unit
                        res[n,c(i)] <-1
                        n<- n+1
                }
                for(j in 4:6){ ## Two units
                        if (num==2){
                                res[n,c(i,j)] <-1
                                n<- n+1
                        }
                        for(k in 7:9){ ## Three units
                                if (num==3){
                                        res[n,c(i,j,k)] <-1
                                        n<- n+1
                                }
                                for(l in 10:12){ ## Four units
                                        if (num==4){
                                                res[n,c(i,j,k,l)] <-1
                                                n<- n+1
                                        }
                                }
                        }
                }
        }
        colnames(res) <- rep(c("Available","Partial","Full_Out"),num)
        res
        
}

## Calculate probabilities for trinomial expansion
## Calls function trinomial_exp to calculate possible states
## Total_MW = total capacity of all units in DUID
## n = number of units in DUID
## prob_out = probability of forced outage
## prob_par = probability of partial outage
## PDR = partial outage derating
## MTTR_O = Mean time to repair for forced outage
## MTTR_P = Mean time to repair for partial outage
## MTTR_D = Mean time to repair default - used for combinations of forced and
##      partial outages to aid in convergence = 0.5 by default
## sig = significance number of decimal places to round probabilities to = 5 by default
## rep_zero = report zero probability states. Default is FALSE


trinom_prob <- function(Total_MW,n, prob_out,prob_par,PDR,MTTR_O,MTTR_P,MTTR_D=0.5,sig=5,rep_zero=FALSE){
        ## Calculate average MW per unit
        Av_MW <- Total_MW /  n 
        ## Calculate average MW per unit when derated
        PDR <- PDR / 100
        P_MW <- Av_MW * (1-PDR)
        
        ## Calculate probability of no outage
        prob_avail <- (1-(prob_out/100))*(1-(prob_par/100))
        ## Calculate probability of a partial outage
        prob_pout <- (1-(prob_out/100))*prob_par/100
       
        prob_out <- (prob_out/100)
        prob <- c(prob_avail, prob_pout,prob_out)
        
        
        ## Calculate range of possibilities by calling the trinomial_exp function
        Comb_State <- trinomial_exp(n)
        ## Assign names to each column in Matrix
        colnames(Comb_State) <- c("Avail","Par","Out")
        
        ## Work out total capcity available in each state (forced outage capacity = 0)
        Comb_Cap <- Comb_State[,1]*Av_MW + Comb_State[,2]*P_MW
        
        ## Probabilities of each operating state
        prob <- c(prob_avail,prob_pout,prob_out)
        
        ## Calculate trinomial probabilties
        Comb_Prob <-(factorial(n)/(factorial(Comb_State[,1])*factorial(Comb_State[,2]) * 
                factorial(Comb_State[,3])))* (prob[1]^Comb_State[,1])*
                (prob[2]^Comb_State[,2])*(prob[3]^Comb_State[,3])
        
        ## Round to level of significance specified in parameters
        Comb_Prob <- round(Comb_Prob*100,5)
        
        ## Calculate Mean Time to Repair. Initialise vector with default values
        MTTR <- rep(MTTR_D,nrow(Comb_State))
        ## Where All units are available except one is on partial outage, use MTTR_P
        MTTR[Comb_State[,"Par"]==1 & Comb_State[,"Avail"]==(n-1)] <- MTTR_P
        ## Where all units are available except for one on full outage, use MTTR_O
        MTTR[Comb_State[,"Out"]==1 & Comb_State[,"Avail"]==(n-1)] <- MTTR_O
        ## Where all units are available, MTTR is zero
        MTTR[Comb_State[,"Avail"]==n] <- 0
        
        ## Calculate Percentage of Capacity on Outage
        CapOut <- round((Total_MW - Comb_Cap)/Total_MW*100,2)
        
        ## Combine all data for reporting
        Comb <- data.frame(cbind(Comb_State,Comb_Cap,CapOut, MTTR, Comb_Prob))
        ## Name columns of Data
        names(Comb) <- c("Available", "Partial", "Full_Out","Total_Cap","Out_Fac","MTTR (hr)","FOR")
        Comb <- Comb[order(Comb$Out_Fac),]
        
        ## Only return probabilities above zero (According to specified significance)
        if (rep_zero==FALSE)
                return(Comb[Comb[,"FOR"]> 0,c("Out_Fac","MTTR (hr)","FOR")])
        else if (rep_zero==TRUE)
                return(Comb)
}



## Full Expansion 
## State = matrix of possible operating states from full_exp
## Prob = vector of probabilities(prob_avail, prob_pout, prob_out)
## unit_cap = vector of individual unit capacities
## n = number of units
## PDR = partial outage derating
## MTTR_P - MTTR Partial Outage MTTR_O = Full Outage, MTTR_D = Default for combinations
## returns a vector with one probability for each state
fullexp_prob <- function(unit_cap,n,prob_out,prob_par,PDR,MTTR_O,MTTR_P,MTTR_D=0.5,CCGT=FALSE,sig=4,rep_zero=FALSE){
        prob_avail <- (1-(prob_out/100))*(1-(prob_par/100))
        prob_pout <- (1-(prob_out/100))*prob_par/100
        prob_out <- (prob_out/100)
        prob <- c(prob_avail, prob_pout,prob_out)
        PDR <- PDR / 100
        
        state <- full_exp(n)
        
        ## count the total number of units in each state
        if(n==1){
                Avail <- state[,1]
                Par <- state[,2]
                Out <- state[,3]
                Comb_Cap <- Avail* unit_cap + (1-PDR)*Par*unit_cap
        }
        
        else if (n==2){
                Avail <- state[,1] + state[,4]
                Avail_Cap <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2]
                Par <- state[,2] + state[,5]
                Par_Cap <- (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2])
                Out <- state[,3] + state[,6]
                Comb_Cap <- Avail_Cap + Par_Cap
                if(CCGT==TRUE){
                        ##No Capacity if Gas Unit is out. Always assume steam unit is last one
                        Comb_Cap[state[,3]==1] <- 0
                }
        }
        else if(n==3){
                Avail <- state[,1] + state[,4] + state[,7]
                Avail_Cap <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2]+ state[,7]*unit_cap[3]
                Par <- state[,2] + state[,5] + state[,8]
                Par_Cap <- (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2]+ state[,8]*unit_cap[3])
                Out <- state[,3] + state[,6] + state[,9]
                
                Comb_Cap <- Avail_Cap + Par_Cap
                if(CCGT==TRUE){
                        ##No Capacity if both Gas Units are out. Half capacity if one out. Always assume steam unit is last one
                        Gas_Out <- state[,3] + state[,6]
                        Gas_On <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2] + (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2])
                        Steam_On <- state[,7]*unit_cap[3]+ (1-PDR)*state[,8]*unit_cap[3]
                        Comb_Cap[(state[,9]==0 & Gas_Out==1)] <- Gas_On[(state[,9]==0 & Gas_Out==1)] + 1/2* Steam_On[(state[,9]==0 & Gas_Out==1)]
                        Comb_Cap[(state[,9]==0 & Gas_Out==2)] <- 0
                }
        }
        else if(n==4){
                Avail <- state[,1] + state[,4] + state[,7] + state[,10]
                Avail_Cap <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2]+ state[,7]*unit_cap[3]+ state[,10]*unit_cap[4]
                Par <- state[,2] + state[,5] + state[,8] + state[,11]
                Par_Cap <- (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2]+ state[,8]*unit_cap[3]+ state[,11]*unit_cap[4])
                Out <- state[,3] + state[,6] + state[,9] + state[,12]
                
                Comb_Cap <- Avail_Cap + Par_Cap
                
                if(CCGT==TRUE){
                        ##No Capacity if both Gas Units are out. Half capacity if one out. Always assume steam unit is last one
                        Gas_Out <- state[,3] + state[,6] + state[,9]
                        Gas_On <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2] + state[,7]*unit_cap[3]+
                                (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2]+ state[,8]*unit_cap[3])
                        Steam_On <- state[,10]*unit_cap[4]+ (1-PDR)*state[,11]*unit_cap[4]
                        Comb_Cap[(state[,12]==0 & Gas_Out==1)] <- Gas_On[(state[,12]==0 & Gas_Out==1)] + 2/3* Steam_On[(state[,12]==0 & Gas_Out==1)]
                        Comb_Cap[(state[,12]==0 & Gas_Out==2)] <- Gas_On[(state[,12]==0 & Gas_Out==2)] + 1/3* Steam_On[(state[,12]==0 & Gas_Out==2)]
                        Comb_Cap[(state[,12]==0 & Gas_Out==3)] <- 0
                }
        }
        Comb_State <- cbind(Avail,Par,Out)
        
        ##Calculate Mean Time To Repair
        
        MTTR <- rep(MTTR_D,3^n)
        MTTR[Comb_State[,"Par"]==1 & Comb_State[,"Avail"]==(n-1)] <- MTTR_P
        MTTR[Comb_State[,"Out"]==1 & Comb_State[,"Avail"]==(n-1)] <- MTTR_O
        MTTR[Comb_State[,"Avail"]==n] <- 0
        
        ## Calculate probabilties
        Comb_Prob <-(prob[1]^Comb_State[,1])*(prob[2]^Comb_State[,2])*(prob[3]^Comb_State[,3])
        Comb_Prob <- round(Comb_Prob*100,sig)
        
        ## Combine all data
        
        ## Calculate Percentage of Capacity on Outage
        CapOut <- round((sum(unit_cap) - Comb_Cap)/sum(unit_cap)*100,2)
        Comb <- data.frame(cbind(Comb_State,Comb_Cap,CapOut, Comb_Prob,MTTR))
        ## Sum up probabilities of equal operating states to reduce overall size
        ## i.e. where the same configuration of operating states has the same available capacity
        Reduced <- aggregate(Comb$Comb_Prob,by=list(Comb$CapOut, Comb$MTTR),FUN=sum)
        names(Reduced) <- c("Out_Fac","MTTR (hr)","FOR")
        
        ## Only return probabilities above zero (5dp)
        Reduced <- Reduced[order(Reduced$Out_Fac),]
        
        ## Only return probabilities above zero (According to specified significance)
        if (rep_zero==FALSE)
                return(Reduced[Reduced[,"FOR"]> 0,])
        else if (rep_zero==TRUE)
                return(Comb)
}


out.rates <- read.csv("outage.csv")
rownames(out.rates) <- as.character(out.rates[,1])
out.rates <- out.rates[,2:6]





Murray <- Outage_Calc(1500,14,0.0082217,0.0000333,0.119697,20.45*2,5.9*2)
Hallet <- Outage_Calc(220,12,0.0066015,0.0008687,0.3189,33.45*2,10.75*2)
SHGEN <- Outage_Calc(240,4,0.0082217,0.0000333,0.119697,20.45*2,5.9*2)

write.to.plexos <- function(GenInfo,DUID,GenType){
        ## Subtract 1 from number of bands as we don't write out Band 1 (Full Capacity)
        num_bands <- nrow(GenInfo)-1
       
        #Initialise output matrix
        plex_out <- matrix(ncol=13,nrow=num_bands*3)
        
        #List of properties to write to Plexos, with fields and units
        properties <- c("Forced Outage Rate","Outage Factor","Mean Time to Repair")
        fields <- c("FOR","Out_Fac","MTTR (hr)")
        units <- c("%","%","h")
        
        #Set row count to 1
        n <- 1
        for (i in 1:3){
                ## Go through one property at a time
                prop <- properties[i]
                unit <- units[i]
                ## Write out all the bands for each property
                for (j in 1:num_bands){
                        field <- fields[i]
                        #Band 1 is ignored as it represents full capacity - start at (j+1)
                        new_row <- c(prop,GenInfo[j+1,field],"",unit,j+1,
                                     "","","","=","","","",GenType)
                        plex_out[n,] <- new_row
                      #Move to next row
                        n <- n+1
                }
        }
        
        ## Insert DUID name as first column
        plex_out <- cbind(DUID,plex_out)
        ## Label output columns
        colnames(plex_out) <- c("DUID","Property","Value","Data File","Units","Band",
                              "Date From","Date To","Timeslice","Action","Expression",
                              "Scenario","Memo","Category")
        data.frame(plex_out)
}


LEM_WIL <- fullexp_prob(unit_cap=c(51,30.6),n=2,
                        prob_out=out.rates["Hydro","P_OUT"],
                        prob_par=out.rates["Hydro","P_PAR"],
                        PDR=out.rates["Hydro","PDR"],
                        MTTR_O=out.rates["Hydro","MTTR_O"],
                        MTTR_P=out.rates["Hydro","MTTR_P"])


SHGEN <- fullexp_prob(unit_cap=c(40,40,80,80),n=4,
                        prob_out=out.rates["Hydro","P_OUT"],
                        prob_par=out.rates["Hydro","P_PAR"],
                        PDR=out.rates["Hydro","PDR"],
                        MTTR_O=out.rates["Hydro","MTTR_O"],
                        MTTR_P=out.rates["Hydro","MTTR_P"])

CPSA <- fullexp_prob(unit_cap=c(43,43,57),n=3,
                      prob_out=out.rates["CCGT","P_OUT"],
                      prob_par=out.rates["CCGT","P_PAR"],
                      PDR=out.rates["CCGT","PDR"],
                      MTTR_O=out.rates["CCGT","MTTR_O"],
                      MTTR_P=out.rates["CCGT","MTTR_P"],
                      CCGT=TRUE)

DDPS1 <- fullexp_prob(unit_cap=c(121,121,121,280),n=4,
                     prob_out=out.rates["CCGT","P_OUT"],
                     prob_par=out.rates["CCGT","P_PAR"],
                     PDR=out.rates["CCGT","PDR"],
                     MTTR_O=out.rates["CCGT","MTTR_O"],
                     MTTR_P=out.rates["CCGT","MTTR_P"],
                     CCGT=TRUE)

OSB_AG <- fullexp_prob(unit_cap=c(118,62),n=2,
                      prob_out=out.rates["CCGT","P_OUT"],
                      prob_par=out.rates["CCGT","P_PAR"],
                      PDR=out.rates["CCGT","PDR"],
                      MTTR_O=out.rates["CCGT","MTTR_O"],
                      MTTR_P=out.rates["CCGT","MTTR_P"],
                      CCGT=TRUE)

PPCCGT <- fullexp_prob(unit_cap=c(160,160,158),n=3,
                     prob_out=out.rates["CCGT","P_OUT"],
                     prob_par=out.rates["CCGT","P_PAR"],
                     PDR=out.rates["CCGT","PDR"],
                     MTTR_O=out.rates["CCGT","MTTR_O"],
                     MTTR_P=out.rates["CCGT","MTTR_P"],
                     CCGT=TRUE)

SITHE01 <- fullexp_prob(unit_cap=c(30,30,30,62),n=4,
                      prob_out=out.rates["CCGT","P_OUT"],
                      prob_par=out.rates["CCGT","P_PAR"],
                      PDR=out.rates["CCGT","PDR"],
                      MTTR_O=out.rates["CCGT","MTTR_O"],
                      MTTR_P=out.rates["CCGT","MTTR_P"],
                      CCGT=TRUE)

TVCC201 <- fullexp_prob(unit_cap=c(68,141),n=2,
                       prob_out=out.rates["CCGT","P_OUT"],
                       prob_par=out.rates["CCGT","P_PAR"],
                       PDR=out.rates["CCGT","PDR"],
                       MTTR_O=out.rates["CCGT","MTTR_O"],
                       MTTR_P=out.rates["CCGT","MTTR_P"],
                       CCGT=TRUE)

AGLHAL <- trinom_prob(Total_MW=222,n=12,
                        prob_out=out.rates["OCGT","P_OUT"],
                        prob_par=out.rates["OCGT","P_PAR"],
                        PDR=out.rates["OCGT","PDR"],
                        MTTR_O=out.rates["OCGT","MTTR_O"],
                        MTTR_P=out.rates["OCGT","MTTR_P"])

AGLSOM <- trinom_prob(Total_MW=160,n=4,
                      prob_out=out.rates["OCGT","P_OUT"],
                      prob_par=out.rates["OCGT","P_PAR"],
                      PDR=out.rates["OCGT","PDR"],
                      MTTR_O=out.rates["OCGT","MTTR_O"],
                      MTTR_P=out.rates["OCGT","MTTR_P"])

HVGTS <- trinom_prob(Total_MW=50,n=2,
                      prob_out=out.rates["OCGT","P_OUT"],
                      prob_par=out.rates["OCGT","P_PAR"],
                      PDR=out.rates["OCGT","PDR"],
                      MTTR_O=out.rates["OCGT","MTTR_O"],
                      MTTR_P=out.rates["OCGT","MTTR_P"])

POR01 <- trinom_prob(Total_MW=50,n=2,
                      prob_out=out.rates["OCGT","P_OUT"],
                      prob_par=out.rates["OCGT","P_PAR"],
                      PDR=out.rates["OCGT","PDR"],
                      MTTR_O=out.rates["OCGT","MTTR_O"],
                      MTTR_P=out.rates["OCGT","MTTR_P"])

SNUG1 <- trinom_prob(Total_MW=63,n=3,
                      prob_out=out.rates["OCGT","P_OUT"],
                      prob_par=out.rates["OCGT","P_PAR"],
                      PDR=out.rates["OCGT","PDR"],
                      MTTR_O=out.rates["OCGT","MTTR_O"],
                      MTTR_P=out.rates["OCGT","MTTR_P"])


GORDON <- trinom_prob(Total_MW=432,n=3,
                      prob_out=out.rates["Hydro","P_OUT"],
                      prob_par=out.rates["Hydro","P_PAR"],
                      PDR=out.rates["Hydro","PDR"],
                      MTTR_O=out.rates["Hydro","MTTR_O"],
                      MTTR_P=out.rates["Hydro","MTTR_P"])


GUTHEGA<- trinom_prob(Total_MW=60,n=2,
                      prob_out=out.rates["Hydro","P_OUT"],
                      prob_par=out.rates["Hydro","P_PAR"],
                      PDR=out.rates["Hydro","PDR"],
                      MTTR_O=out.rates["Hydro","MTTR_O"],
                      MTTR_P=out.rates["Hydro","MTTR_P"])

LI_WY_CA<- trinom_prob(Total_MW=183,n=8,
                       prob_out=out.rates["Hydro","P_OUT"],
                       prob_par=out.rates["Hydro","P_PAR"],
                       PDR=out.rates["Hydro","PDR"],
                       MTTR_O=out.rates["Hydro","MTTR_O"],
                       MTTR_P=out.rates["Hydro","MTTR_P"])

MCKAY1<- trinom_prob(Total_MW=300,n=8,
                     prob_out=out.rates["Hydro","P_OUT"],
                     prob_par=out.rates["Hydro","P_PAR"],
                     PDR=out.rates["Hydro","PDR"],
                     MTTR_O=out.rates["Hydro","MTTR_O"],
                     MTTR_P=out.rates["Hydro","MTTR_P"])

MURRAY<- trinom_prob(Total_MW=1500,n=14,
                     prob_out=out.rates["Hydro","P_OUT"],
                     prob_par=out.rates["Hydro","P_PAR"],
                     PDR=out.rates["Hydro","PDR"],
                     MTTR_O=out.rates["Hydro","MTTR_O"],
                     MTTR_P=out.rates["Hydro","MTTR_P"])

POAT110<- trinom_prob(Total_MW=100,n=2,
                      prob_out=out.rates["Hydro","P_OUT"],
                      prob_par=out.rates["Hydro","P_PAR"],
                      PDR=out.rates["Hydro","PDR"],
                      MTTR_O=out.rates["Hydro","MTTR_O"],
                      MTTR_P=out.rates["Hydro","MTTR_P"])

POAT220<- trinom_prob(Total_MW=200,n=4,
                      prob_out=out.rates["Hydro","P_OUT"],
                      prob_par=out.rates["Hydro","P_PAR"],
                      PDR=out.rates["Hydro","PDR"],
                      MTTR_O=out.rates["Hydro","MTTR_O"],
                      MTTR_P=out.rates["Hydro","MTTR_P"])

TARRALEA<- trinom_prob(Total_MW=90,n=6,
                       prob_out=out.rates["Hydro","P_OUT"],
                       prob_par=out.rates["Hydro","P_PAR"],
                       PDR=out.rates["Hydro","PDR"],
                       MTTR_O=out.rates["Hydro","MTTR_O"],
                       MTTR_P=out.rates["Hydro","MTTR_P"])

TREVALLN<- trinom_prob(Total_MW=80,n=4,
                       prob_out=out.rates["Hydro","P_OUT"],
                       prob_par=out.rates["Hydro","P_PAR"],
                       PDR=out.rates["Hydro","PDR"],
                       MTTR_O=out.rates["Hydro","MTTR_O"],
                       MTTR_P=out.rates["Hydro","MTTR_P"])

TUMUT3<- trinom_prob(Total_MW=1500,n=6,
                     prob_out=out.rates["Hydro","P_OUT"],
                     prob_par=out.rates["Hydro","P_PAR"],
                     PDR=out.rates["Hydro","PDR"],
                     MTTR_O=out.rates["Hydro","MTTR_O"],
                     MTTR_P=out.rates["Hydro","MTTR_P"])

TUNGATIN<- trinom_prob(Total_MW=125,n=5,
                       prob_out=out.rates["Hydro","P_OUT"],
                       prob_par=out.rates["Hydro","P_PAR"],
                       PDR=out.rates["Hydro","PDR"],
                       MTTR_O=out.rates["Hydro","MTTR_O"],
                       MTTR_P=out.rates["Hydro","MTTR_P"])

UPPTUMUT<- trinom_prob(Total_MW=616,n=8,
                       prob_out=out.rates["Hydro","P_OUT"],
                       prob_par=out.rates["Hydro","P_PAR"],
                       PDR=out.rates["Hydro","PDR"],
                       MTTR_O=out.rates["Hydro","MTTR_O"],
                       MTTR_P=out.rates["Hydro","MTTR_P"])


ANGAST1<- trinom_prob(Total_MW=50,n=30,
                       prob_out=out.rates["Gas_other","P_OUT"],
                       prob_par=out.rates["Gas_other","P_PAR"],
                       PDR=out.rates["Gas_other","PDR"],
                       MTTR_O=out.rates["Gas_other","MTTR_O"],
                       MTTR_P=out.rates["Gas_other","MTTR_P"])

LONSDALE<- trinom_prob(Total_MW=20,n=18,
                      prob_out=out.rates["Gas_other","P_OUT"],
                      prob_par=out.rates["Gas_other","P_PAR"],
                      PDR=out.rates["Gas_other","PDR"],
                      MTTR_O=out.rates["Gas_other","MTTR_O"],
                      MTTR_P=out.rates["Gas_other","MTTR_P"])

fullset <- rbind(write.to.plexos(LEM_WIL,"LEM_WIL","Hydro"),
                 write.to.plexos(SHGEN,"SHGEN","Hydro"),
                 write.to.plexos(CPSA,"CPSA","CCGT"),
                 write.to.plexos(DDPS1,"DDPS1","CCGT"),
                 write.to.plexos(OSB_AG,"OSB-AG","CCGT"),
                 write.to.plexos(PPCCGT,"PPCCGT","CCGT"),
                 write.to.plexos(SITHE01,"SITHE01","CCGT"),
                 write.to.plexos(TVCC201,"TVCC201","CCGT"),
                 write.to.plexos(ANGAST1,"ANGAST1","Diesel"),
                 write.to.plexos(LONSDALE,"LONSDALE","Diesel"),
                 write.to.plexos(GORDON,"GORDON","Hydro"),
                 write.to.plexos(GUTHEGA,"GUTHEGA","Hydro"),
                 write.to.plexos(LI_WY_CA,"LI_WY_CA","Hydro"),
                 write.to.plexos(MCKAY1,"MCKAY1","Hydro"),
                 write.to.plexos(MURRAY,"MURRAY","Hydro"),
                 write.to.plexos(POAT110,"POAT110","Hydro"),
                 write.to.plexos(POAT220,"POAT220","Hydro"),
                 write.to.plexos(TARRALEA,"TARRALEA","Hydro"),
                 write.to.plexos(TREVALLN,"TREVALLN","Hydro"),
                 write.to.plexos(TUMUT3,"TUMUT3","Hydro"),
                 write.to.plexos(TUNGATIN,"TUNGATIN","Hydro"),
                 write.to.plexos(UPPTUMUT,"UPPTUMUT","Hydro"),
                 write.to.plexos(AGLHAL,"AGLHAL","OCGT"),
                 write.to.plexos(AGLSOM,"AGLSOM","OCGT"),
                 write.to.plexos(HVGTS,"HVGTS","OCGT"),
                 write.to.plexos(POR01,"POR01","OCGT"),
                 write.to.plexos(SNUG1,"SNUG1","OCGT"))
                 
fullset[,"Value"] <- as.numeric(as.character(fullset[,"Value"]))
