## Read in outage rates

out.rates <- read.csv("outage.csv")
rownames(out.rates) <- as.character(out.rates[,1])
out.rates <- out.rates[,2:6]

## FUNCTION
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

## FUNCTION
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

## FUNCTION trinom_prob
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


trinom_prob <- function(Total_MW,n, out.rates,type, MTTR_D=0.5,sig=5,rep_zero=FALSE){
        
        # Extract probabililties and MTTR from out.rates
        prob_out <- out.rates[type,"P_OUT"]
        prob_par <- out.rates[type,"P_PAR"]
        PDR <- out.rates[type,"PDR"]/100
        MTTR_O <- out.rates[type,"MTTR_O"]
        MTTR_P <- out.rates[type,"MTTR_P"]
        ## Calculate average MW per unit
        Av_MW <- Total_MW /  n
        ## Calculate average MW per unit when derated
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
        #return(Comb[Comb[,"FOR"]> 0,c("Available", "Partial", "Full_Out","Total_Cap","FOR")])
        else if (rep_zero==TRUE)
                return(Comb)
}


## FUNCTION fullexp_prob
## Full Expansion of all possible states for n = 1, 2, 3 or 4 only
## unit_cap = vector of individual unit capacities
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
fullexp_prob <- function(unit_cap,n, out.rates,type, MTTR_D=0.5,CCGT=FALSE,sig=5,rep_zero=FALSE){
        
        ## Calculate the probabilities of each state
        ## They are entered as values between 0 and 100, so need to divide by 100
        
        # Extract probabililties and MTTR from out.rates
        prob_out <- out.rates[type,"P_OUT"]/100
        prob_par <- out.rates[type,"P_PAR"]/100
        PDR <- out.rates[type,"PDR"]/100
        MTTR_O <- out.rates[type,"MTTR_O"]
        MTTR_P <- out.rates[type,"MTTR_P"]
        
        ## Probability that a unit is available
        prob_avail <- (1-(prob_out))*(1-(prob_par))
        ## Probability that a unit is partially out
        prob_pout <- (1-(prob_out))*prob_par
        
        
        # Vector containing probabilities
        prob <- c(prob_avail, prob_pout,prob_out)
        
        ## Call function to generate matrix of all possible states
        state <- full_exp(n)
        
        ## Different calculations depending on number of units in DUID
        ## If DUID is a CCGT and n >=2, assume last unit is the steam unit
        ## and calculate operating states for the
        ## Single Unit
        if(n==1){
                Avail <- state[,1]
                Par <- state[,2]
                Out <- state[,3]
                ## Calculate Capacity of DUID for each state.
                ## No capacity if if on forced outage.
                Comb_Cap <- Avail* unit_cap + (1-PDR)*Par*unit_cap
        }
        ## Two Units
        else if (n==2){
                ## Number of units available
                Avail <- state[,1] + state[,4]
                ## Combined capacity of available units
                Avail_Cap <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2]
                ## Number of units partially out
                Par <- state[,2] + state[,5]
                ## Capacity of units partially out
                Par_Cap <- (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2])
                ## Number of units on full outage
                Out <- state[,3] + state[,6]
                ## Calculate total capacity of DUID
                Comb_Cap <- Avail_Cap + Par_Cap
                
                ## Calculate values for CCGT
                if(CCGT==TRUE){
                        ## Combined Capacity of DUID is 0 if Gas Unit is out
                        Comb_Cap[state[,3]==1] <- 0
                }
        }
        ## Three Units
        else if(n==3){
                ##Number of units available
                Avail <- state[,1] + state[,4] + state[,7]
                ## Combined capacity of available units
                Avail_Cap <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2]+
                        state[,7]*unit_cap[3]
                ## Number of units partially out
                Par <- state[,2] + state[,5] + state[,8]
                ## Capacity of units partially out
                Par_Cap <- (1-PDR)*(state[,2]*unit_cap[1] +
                                            state[,5]*unit_cap[2]+ state[,8]*unit_cap[3])
                ## Number of units on full outage
                Out <- state[,3] + state[,6] + state[,9]
                ## Calculate total capacity of DUID
                Comb_Cap <- Avail_Cap + Par_Cap
                
                ## Calculate values for CCGT
                if(CCGT==TRUE){
                        ## Calculate how many gas units are out
                        Gas_Out <- state[,3] + state[,6]
                        ## Calculate capacity of gas units available or
                        ## partially operating
                        Gas_On <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2] +
                                (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2])
                        ## Calculate capacity of steam unit
                        Steam_On <- state[,7]*unit_cap[3]+
                                (1-PDR)*state[,8]*unit_cap[3]
                        ## If steam available or partially available with
                        ## one gas out, the capacity is equal to
                        ## gas capacity + 1/2 steam sapacity
                        Comb_Cap[(state[,9]==0 & Gas_Out==1)] <-
                                Gas_On[(state[,9]==0 & Gas_Out==1)] +
                                1/2* Steam_On[(state[,9]==0 & Gas_Out==1)]
                        ## If Steam available and both gas out, capacity is 0
                        Comb_Cap[(state[,9]==0 & Gas_Out==2)] <- 0
                }
        }
        ## Four Units
        else if(n==4){
                ##Number of units available
                Avail <- state[,1] + state[,4] + state[,7] + state[,10]
                ## Combined capacity of available units
                Avail_Cap <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2]+
                        state[,7]*unit_cap[3]+ state[,10]*unit_cap[4]
                ## Number of units partially out
                Par <- state[,2] + state[,5] + state[,8] + state[,11]
                ## Capacity of units partially out
                Par_Cap <- (1-PDR)*(state[,2]*unit_cap[1] + state[,5]*unit_cap[2]+
                                            state[,8]*unit_cap[3]+ state[,11]*unit_cap[4])
                ## Number of units on full outage
                Out <- state[,3] + state[,6] + state[,9] + state[,12]
                ## Calculate total capacity of DUID
                Comb_Cap <- Avail_Cap + Par_Cap
                
                ## Calculate values for CCGT
                if(CCGT==TRUE){
                        ## Calculate how many gas units are out
                        Gas_Out <- state[,3] + state[,6] + state[,9]
                        ## Calculate capacity of gas units available or
                        ## partially operating
                        Gas_On <- state[,1]*unit_cap[1] + state[,4]*unit_cap[2] +
                                state[,7]*unit_cap[3] +
                                (1-PDR)*(state[,2]*unit_cap[1] +
                                                 state[,5]*unit_cap[2]+ state[,8]*unit_cap[3])
                        ## Calculate capacity of steam unit
                        Steam_On <- state[,10]*unit_cap[4]+
                                (1-PDR)*state[,11]*unit_cap[4]
                        ## If steam available or partially available with
                        ## one gas out, the capacity is equal to
                        ## gas capacity + 2/3 steam sapacity
                        Comb_Cap[(state[,12]==0 & Gas_Out==1)] <-
                                Gas_On[(state[,12]==0 & Gas_Out==1)] +
                                2/3* Steam_On[(state[,12]==0 & Gas_Out==1)]
                        ## If steam available or partially available with
                        ## one gas out, the capacity is equal to
                        ## gas capacity + 1/3 steam sapacity
                        Comb_Cap[(state[,12]==0 & Gas_Out==2)] <-
                                Gas_On[(state[,12]==0 & Gas_Out==2)] +
                                1/3* Steam_On[(state[,12]==0 & Gas_Out==2)]
                        ## If Steam available and three gas out, capacity is 0
                        Comb_Cap[(state[,12]==0 & Gas_Out==3)] <- 0
                }
        }
        ## Combined Availability, Partial Outage and Full Outage results
        Comb_State <- cbind(Avail,Par,Out)
        
        ## Calculate Mean Time to Repair. Initialise vector with default values
        MTTR <- rep(MTTR_D,3^n)
        ## Where All units are available except one is on partial outage, use MTTR_P
        MTTR[Comb_State[,"Par"]==1 & Comb_State[,"Avail"]==(n-1)] <- MTTR_P
        ## Where All units are available except one is on full outage, use MTTR_O
        MTTR[Comb_State[,"Out"]==1 & Comb_State[,"Avail"]==(n-1)] <- MTTR_O
        ## Where all units are available, MTTR is zero
        MTTR[Comb_State[,"Avail"]==n] <- 0
        
        ## Calculate probabilties
        Comb_Prob <-(prob[1]^Comb_State[,1])*(prob[2]^Comb_State[,2])*(prob[3]^Comb_State[,3])
        ## Round to level of significance specified
        Comb_Prob <- round(Comb_Prob*100,sig)
        
        ## Calculate Percentage of Capacity on Outage
        CapOut <- round((sum(unit_cap) - Comb_Cap)/sum(unit_cap)*100,2)
        
        ## Combine all results
        Comb <- data.frame(cbind(Comb_State,Comb_Cap,CapOut, Comb_Prob,MTTR))
        
        ## Sum up probabilities of equal operating states to reduce overall size
        ## i.e. where the same configuration of operating states has the same available capacity
        Reduced <- aggregate(Comb$Comb_Prob,by=list(Comb$CapOut, Comb$MTTR),FUN=sum)
        #Reduced <- aggregate(Comb$Comb_Prob,by=list(Comb$Avail,Comb$Par, Comb$Out, Comb$Comb_Cap,Comb$CapOut, Comb$MTTR),FUN=sum)
        ## Label columns
        names(Reduced) <- c("Out_Fac","MTTR (hr)","FOR")
        
        ## Sort results by Outage Factor
        Reduced <- Reduced[order(Reduced$Out_Fac),]
        
        ## Only return probabilities above zero (According to specified significance)
        if (rep_zero==FALSE)
                return(Reduced[Reduced[,"FOR"]> 0,])
        
        ## Report full results
        else if (rep_zero==TRUE)
                return(Comb)
}



write.to.plexos <- function(GenInfo,DUID){
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
                                     "","","","=","","","","")
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

CPSA <- fullexp_prob(unit_cap=c(43,43,57),n=3,out.rates,"CCGT", CCGT=TRUE)
DDPS1 <- fullexp_prob(unit_cap=c(121,121,121,280),n=4,out.rates,"CCGT", CCGT=TRUE)
LEM_WIL <- fullexp_prob(unit_cap=c(51,30.6),n=2,out.rates,"Hydro")
OSB_AG <- fullexp_prob(unit_cap=c(118,62),n=2,out.rates,"CCGT", CCGT=TRUE)
PPCCGT <- fullexp_prob(unit_cap=c(160,160,158),n=3,out.rates,"CCGT", CCGT=TRUE)
SHGEN <- fullexp_prob(unit_cap=c(40,40,80,80),n=4,out.rates,"Hydro")
SITHE01 <- fullexp_prob(unit_cap=c(30,30,30,62),n=4,out.rates,"CCGT", CCGT=TRUE)
TVCC201 <- fullexp_prob(unit_cap=c(68,141),n=2,out.rates,"CCGT", CCGT=TRUE)
GSTONE1 <- trinom_prob(Total_MW=280, n=1,out.rates,"Black_Coal")
GSTONE2 <- trinom_prob(Total_MW=280, n=1,out.rates,"Black_Coal")
GSTONE3 <- trinom_prob(Total_MW=280, n=1,out.rates,"Black_Coal")
GSTONE4 <- trinom_prob(Total_MW=280, n=1,out.rates,"Black_Coal")
GSTONE5 <- trinom_prob(Total_MW=280, n=1,out.rates,"Black_Coal")
GSTONE6 <- trinom_prob(Total_MW=280, n=1,out.rates,"Black_Coal")
CALL_B_1 <- trinom_prob(Total_MW=350, n=1,out.rates,"Black_Coal")
CALL_B_2 <- trinom_prob(Total_MW=350, n=1,out.rates,"Black_Coal")
TARONG1 <- trinom_prob(Total_MW=350, n=1,out.rates,"Black_Coal")
TARONG2 <- trinom_prob(Total_MW=350, n=1,out.rates,"Black_Coal")
TARONG3 <- trinom_prob(Total_MW=350, n=1,out.rates,"Black_Coal")
TARONG4 <- trinom_prob(Total_MW=350, n=1,out.rates,"Black_Coal")
STAN1 <- trinom_prob(Total_MW=365, n=1,out.rates,"Black_Coal")
STAN2 <- trinom_prob(Total_MW=365, n=1,out.rates,"Black_Coal")
STAN3 <- trinom_prob(Total_MW=365, n=1,out.rates,"Black_Coal")
STAN4 <- trinom_prob(Total_MW=365, n=1,out.rates,"Black_Coal")
CPP_3 <- trinom_prob(Total_MW=420, n=1,out.rates,"Black_Coal")
CPP_4 <- trinom_prob(Total_MW=420, n=1,out.rates,"Black_Coal")
MPP_1 <- trinom_prob(Total_MW=426, n=1,out.rates,"Black_Coal")
MPP_2 <- trinom_prob(Total_MW=426, n=1,out.rates,"Black_Coal")
TNPS1 <- trinom_prob(Total_MW=450, n=1,out.rates,"Black_Coal")
BW01 <- trinom_prob(Total_MW=660, n=1,out.rates,"Black_Coal")
BW02 <- trinom_prob(Total_MW=660, n=1,out.rates,"Black_Coal")
BW03 <- trinom_prob(Total_MW=660, n=1,out.rates,"Black_Coal")
BW04 <- trinom_prob(Total_MW=660, n=1,out.rates,"Black_Coal")
VP5 <- trinom_prob(Total_MW=660, n=1,out.rates,"Black_Coal")
VP6 <- trinom_prob(Total_MW=660, n=1,out.rates,"Black_Coal")
MP1 <- trinom_prob(Total_MW=700, n=1,out.rates,"Black_Coal")
MP2 <- trinom_prob(Total_MW=700, n=1,out.rates,"Black_Coal")
ER01 <- trinom_prob(Total_MW=720, n=1,out.rates,"Black_Coal")
ER02 <- trinom_prob(Total_MW=720, n=1,out.rates,"Black_Coal")
ER03 <- trinom_prob(Total_MW=720, n=1,out.rates,"Black_Coal")
ER04 <- trinom_prob(Total_MW=720, n=1,out.rates,"Black_Coal")
KPP_1 <- trinom_prob(Total_MW=744, n=1,out.rates,"Black_Coal")
YWPS1 <- trinom_prob(Total_MW=360, n=1,out.rates,"Brown_Coal")
YWPS2 <- trinom_prob(Total_MW=360, n=1,out.rates,"Brown_Coal")
YWPS3 <- trinom_prob(Total_MW=380, n=1,out.rates,"Brown_Coal")
YWPS4 <- trinom_prob(Total_MW=380, n=1,out.rates,"Brown_Coal")
LOYYB1 <- trinom_prob(Total_MW=500, n=1,out.rates,"Brown_Coal")
LOYYB2 <- trinom_prob(Total_MW=500, n=1,out.rates,"Brown_Coal")
LYA2 <- trinom_prob(Total_MW=530, n=1,out.rates,"Brown_Coal")
LYA1 <- trinom_prob(Total_MW=560, n=1,out.rates,"Brown_Coal")
LYA3 <- trinom_prob(Total_MW=560, n=1,out.rates,"Brown_Coal")
LYA4 <- trinom_prob(Total_MW=560, n=1,out.rates,"Brown_Coal")
BARCALDN <- trinom_prob(Total_MW=37, n=1,out.rates,"CCGT")
SWAN_E <- trinom_prob(Total_MW=385, n=1,out.rates,"CCGT")
ROMA_7 <- trinom_prob(Total_MW=40, n=1,out.rates,"OCGT")
TORRA1 <- trinom_prob(Total_MW=120, n=1,out.rates,"Gas_other")
TORRA2 <- trinom_prob(Total_MW=120, n=1,out.rates,"Gas_other")
TORRA3 <- trinom_prob(Total_MW=120, n=1,out.rates,"Gas_other")
TORRA4 <- trinom_prob(Total_MW=120, n=1,out.rates,"Gas_other")
TORRB1 <- trinom_prob(Total_MW=200, n=1,out.rates,"Gas_other")
TORRB2 <- trinom_prob(Total_MW=200, n=1,out.rates,"Gas_other")
TORRB3 <- trinom_prob(Total_MW=200, n=1,out.rates,"Gas_other")
TORRB4 <- trinom_prob(Total_MW=200, n=1,out.rates,"Gas_other")
NPS <- trinom_prob(Total_MW=500, n=1,out.rates,"Gas_other")
WKIEWA1 <- trinom_prob(Total_MW=15, n=1,out.rates,"Hydro")
WKIEWA2 <- trinom_prob(Total_MW=15, n=1,out.rates,"Hydro")
KAREEYA1 <- trinom_prob(Total_MW=21, n=1,out.rates,"Hydro")
KAREEYA2 <- trinom_prob(Total_MW=21, n=1,out.rates,"Hydro")
KAREEYA3 <- trinom_prob(Total_MW=21, n=1,out.rates,"Hydro")
KAREEYA4 <- trinom_prob(Total_MW=21, n=1,out.rates,"Hydro")
HUMENSW <- trinom_prob(Total_MW=29, n=1,out.rates,"Hydro")
HUMEV <- trinom_prob(Total_MW=29, n=1,out.rates,"Hydro")
BARRON1 <- trinom_prob(Total_MW=30, n=1,out.rates,"Hydro")
BARRON2 <- trinom_prob(Total_MW=30, n=1,out.rates,"Hydro")
LK_ECHO <- trinom_prob(Total_MW=32.4, n=1,out.rates,"Hydro")
MEADOWBK <- trinom_prob(Total_MW=40, n=1,out.rates,"Hydro")
FISHER <- trinom_prob(Total_MW=43.2, n=1,out.rates,"Hydro")
EILDON1 <- trinom_prob(Total_MW=60, n=1,out.rates,"Hydro")
EILDON2 <- trinom_prob(Total_MW=60, n=1,out.rates,"Hydro")
DEVILS_G <- trinom_prob(Total_MW=60, n=1,out.rates,"Hydro")
BLOWERNG <- trinom_prob(Total_MW=70, n=1,out.rates,"Hydro")
BASTYAN <- trinom_prob(Total_MW=79.9, n=1,out.rates,"Hydro")
MACKNTSH <- trinom_prob(Total_MW=79.9, n=1,out.rates,"Hydro")
TRIBUTE <- trinom_prob(Total_MW=82.8, n=1,out.rates,"Hydro")
CETHANA <- trinom_prob(Total_MW=85, n=1,out.rates,"Hydro")
REECE1 <- trinom_prob(Total_MW=115.6, n=1,out.rates,"Hydro")
REECE2 <- trinom_prob(Total_MW=115.6, n=1,out.rates,"Hydro")
JBUTTERS <- trinom_prob(Total_MW=144, n=1,out.rates,"Hydro")
DARTM1 <- trinom_prob(Total_MW=150, n=1,out.rates,"Hydro")
LD01 <- trinom_prob(Total_MW=500, n=1,out.rates,"Liddell")
LD02 <- trinom_prob(Total_MW=500, n=1,out.rates,"Liddell")
LD03 <- trinom_prob(Total_MW=500, n=1,out.rates,"Liddell")
LD04 <- trinom_prob(Total_MW=500, n=1,out.rates,"Liddell")
PTSTAN1 <- trinom_prob(Total_MW=1.6, n=1,out.rates,"OCGT")
POR03 <- trinom_prob(Total_MW=23.5, n=1,out.rates,"OCGT")
QPS1 <- trinom_prob(Total_MW=24, n=1,out.rates,"OCGT")
QPS2 <- trinom_prob(Total_MW=24, n=1,out.rates,"OCGT")
QPS3 <- trinom_prob(Total_MW=24, n=1,out.rates,"OCGT")
QPS4 <- trinom_prob(Total_MW=24, n=1,out.rates,"OCGT")
MACKAYGT <- trinom_prob(Total_MW=30, n=1,out.rates,"OCGT")
BBTHREE1 <- trinom_prob(Total_MW=35, n=1,out.rates,"OCGT")
BBTHREE2 <- trinom_prob(Total_MW=35, n=1,out.rates,"OCGT")
BBTHREE3 <- trinom_prob(Total_MW=35, n=1,out.rates,"OCGT")
LADBROK1 <- trinom_prob(Total_MW=40, n=1,out.rates,"OCGT")
LADBROK2 <- trinom_prob(Total_MW=40, n=1,out.rates,"OCGT")
ROMA_8 <- trinom_prob(Total_MW=40, n=1,out.rates,"OCGT")
BDL01 <- trinom_prob(Total_MW=47, n=1,out.rates,"OCGT")
BDL02 <- trinom_prob(Total_MW=47, n=1,out.rates,"OCGT")
VPGS1 <- trinom_prob(Total_MW=50, n=1,out.rates,"OCGT")
VPGS2 <- trinom_prob(Total_MW=50, n=1,out.rates,"OCGT")
VPGS3 <- trinom_prob(Total_MW=50, n=1,out.rates,"OCGT")
VPGS4 <- trinom_prob(Total_MW=50, n=1,out.rates,"OCGT")
VPGS5 <- trinom_prob(Total_MW=50, n=1,out.rates,"OCGT")
VPGS6 <- trinom_prob(Total_MW=50, n=1,out.rates,"OCGT")
JLA01 <- trinom_prob(Total_MW=51, n=1,out.rates,"OCGT")
JLA02 <- trinom_prob(Total_MW=51, n=1,out.rates,"OCGT")
JLA03 <- trinom_prob(Total_MW=51, n=1,out.rates,"OCGT")
JLA04 <- trinom_prob(Total_MW=51, n=1,out.rates,"OCGT")
DRYCGT1 <- trinom_prob(Total_MW=52, n=1,out.rates,"OCGT")
DRYCGT2 <- trinom_prob(Total_MW=52, n=1,out.rates,"OCGT")
DRYCGT3 <- trinom_prob(Total_MW=52, n=1,out.rates,"OCGT")
TVPP104 <- trinom_prob(Total_MW=58, n=1,out.rates,"OCGT")
JLB01 <- trinom_prob(Total_MW=76, n=1,out.rates,"OCGT")
JLB02 <- trinom_prob(Total_MW=76, n=1,out.rates,"OCGT")
JLB03 <- trinom_prob(Total_MW=76, n=1,out.rates,"OCGT")
YABULU2 <- trinom_prob(Total_MW=82, n=1,out.rates,"OCGT")
MINTARO <- trinom_prob(Total_MW=90, n=1,out.rates,"OCGT")
QPS5 <- trinom_prob(Total_MW=128, n=1,out.rates,"OCGT")
MSTUART3 <- trinom_prob(Total_MW=131, n=1,out.rates,"OCGT")
OAKEY1 <- trinom_prob(Total_MW=141, n=1,out.rates,"OCGT")
OAKEY2 <- trinom_prob(Total_MW=141, n=1,out.rates,"OCGT")
MSTUART1 <- trinom_prob(Total_MW=146, n=1,out.rates,"OCGT")
MSTUART2 <- trinom_prob(Total_MW=146, n=1,out.rates,"OCGT")
LNGS1 <- trinom_prob(Total_MW=156, n=1,out.rates,"OCGT")
LNGS2 <- trinom_prob(Total_MW=156, n=1,out.rates,"OCGT")
YABULU <- trinom_prob(Total_MW=160, n=1,out.rates,"OCGT")
URANQ11 <- trinom_prob(Total_MW=166, n=1,out.rates,"OCGT")
URANQ12 <- trinom_prob(Total_MW=166, n=1,out.rates,"OCGT")
URANQ13 <- trinom_prob(Total_MW=166, n=1,out.rates,"OCGT")
URANQ14 <- trinom_prob(Total_MW=166, n=1,out.rates,"OCGT")
BRAEMAR1 <- trinom_prob(Total_MW=168, n=1,out.rates,"OCGT")
BRAEMAR2 <- trinom_prob(Total_MW=168, n=1,out.rates,"OCGT")
BRAEMAR3 <- trinom_prob(Total_MW=168, n=1,out.rates,"OCGT")
BRAEMAR5 <- trinom_prob(Total_MW=173, n=1,out.rates,"OCGT")
BRAEMAR6 <- trinom_prob(Total_MW=173, n=1,out.rates,"OCGT")
BRAEMAR7 <- trinom_prob(Total_MW=173, n=1,out.rates,"OCGT")
CG1 <- trinom_prob(Total_MW=181, n=1,out.rates,"OCGT")
CG2 <- trinom_prob(Total_MW=181, n=1,out.rates,"OCGT")
CG3 <- trinom_prob(Total_MW=181, n=1,out.rates,"OCGT")
CG4 <- trinom_prob(Total_MW=181, n=1,out.rates,"OCGT")
MORTLK11 <- trinom_prob(Total_MW=283, n=1,out.rates,"OCGT")
MORTLK12 <- trinom_prob(Total_MW=283, n=1,out.rates,"OCGT")
TALWA1 <- trinom_prob(Total_MW=460, n=1,out.rates,"OCGT")
AGLHAL <- trinom_prob(Total_MW=222,n=12,out.rates,"OCGT")
AGLSOM <- trinom_prob(Total_MW=160,n=4,out.rates,"OCGT")
ANGAST1<- trinom_prob(Total_MW=50,n=2,out.rates,"OCGT")
GORDON <- trinom_prob(Total_MW=432,n=3,out.rates,"Hydro")
GUTHEGA<- trinom_prob(Total_MW=60,n=2,out.rates,"Hydro")
HVGTS <- trinom_prob(Total_MW=50,n=2,out.rates,"OCGT")
LI_WY_CA<- trinom_prob(Total_MW=183,n=8,out.rates,"Hydro")
LONSDALE<- trinom_prob(Total_MW=20,n=1,out.rates,"OCGT")
MCKAY1<- trinom_prob(Total_MW=300,n=8,out.rates,"Hydro")
MURRAY<- trinom_prob(Total_MW=1500,n=14,out.rates,"Hydro")
POAT110<- trinom_prob(Total_MW=100,n=2,out.rates,"Hydro")
POAT220<- trinom_prob(Total_MW=200,n=4,out.rates,"Hydro")
POR01 <- trinom_prob(Total_MW=50,n=2,out.rates,"OCGT")
SNUG1 <- trinom_prob(Total_MW=63,n=3,out.rates,"OCGT")
TARRALEA<- trinom_prob(Total_MW=90,n=6,out.rates,"Hydro")
TREVALLN<- trinom_prob(Total_MW=80,n=4,out.rates,"Hydro")
TUMUT3<- trinom_prob(Total_MW=1500,n=6,out.rates,"Hydro")
TUNGATIN<- trinom_prob(Total_MW=125,n=5,out.rates,"Hydro")
UPPTUMUT<- trinom_prob(Total_MW=616,n=8,out.rates,"Hydro")



fullset <- rbind(write.to.plexos(CPSA, "CPSA"),
                 write.to.plexos(DDPS1, "DDPS1"),
                 write.to.plexos(LEM_WIL, "LEM_WIL"),
                 write.to.plexos(OSB_AG, "OSB-AG"),
                 write.to.plexos(PPCCGT, "PPCCGT"),
                 write.to.plexos(SHGEN , "SHGEN "),
                 write.to.plexos(SITHE01, "SITHE01"),
                 write.to.plexos(TVCC201, "TVCC201"),
                 write.to.plexos(GSTONE1, "GSTONE1"),
                 write.to.plexos(GSTONE2, "GSTONE2"),
                 write.to.plexos(GSTONE3, "GSTONE3"),
                 write.to.plexos(GSTONE4, "GSTONE4"),
                 write.to.plexos(GSTONE5, "GSTONE5"),
                 write.to.plexos(GSTONE6, "GSTONE6"),
                 write.to.plexos(CALL_B_1, "CALL_B_1"),
                 write.to.plexos(CALL_B_2, "CALL_B_2"),
                 write.to.plexos(TARONG1, "TARONG#1"),
                 write.to.plexos(TARONG2, "TARONG#2"),
                 write.to.plexos(TARONG3, "TARONG#3"),
                 write.to.plexos(TARONG4, "TARONG#4"),
                 write.to.plexos(STAN1, "STAN-1"),
                 write.to.plexos(STAN2, "STAN-2"),
                 write.to.plexos(STAN3, "STAN-3"),
                 write.to.plexos(STAN4, "STAN-4"),
                 write.to.plexos(CPP_3, "CPP_3"),
                 write.to.plexos(CPP_4, "CPP_4"),
                 write.to.plexos(MPP_1, "MPP_1"),
                 write.to.plexos(MPP_2, "MPP_2"),
                 write.to.plexos(TNPS1, "TNPS1"),
                 write.to.plexos(BW01, "BW01"),
                 write.to.plexos(BW02, "BW02"),
                 write.to.plexos(BW03, "BW03"),
                 write.to.plexos(BW04, "BW04"),
                 write.to.plexos(VP5, "VP5"),
                 write.to.plexos(VP6, "VP6"),
                 write.to.plexos(MP1, "MP1"),
                 write.to.plexos(MP2, "MP2"),
                 write.to.plexos(ER01, "ER01"),
                 write.to.plexos(ER02, "ER02"),
                 write.to.plexos(ER03, "ER03"),
                 write.to.plexos(ER04, "ER04"),
                 write.to.plexos(KPP_1, "KPP_1"),
                 write.to.plexos(YWPS1, "YWPS1"),
                 write.to.plexos(YWPS2, "YWPS2"),
                 write.to.plexos(YWPS3, "YWPS3"),
                 write.to.plexos(YWPS4, "YWPS4"),
                 write.to.plexos(LOYYB1, "LOYYB1"),
                 write.to.plexos(LOYYB2, "LOYYB2"),
                 write.to.plexos(LYA2, "LYA2"),
                 write.to.plexos(LYA1, "LYA1"),
                 write.to.plexos(LYA3, "LYA3"),
                 write.to.plexos(LYA4, "LYA4"),
                 write.to.plexos(BARCALDN, "BARCALDN"),
                 write.to.plexos(SWAN_E, "SWAN_E"),
                 write.to.plexos(ROMA_7, "ROMA_7"),
                 write.to.plexos(TORRA1, "TORRA1"),
                 write.to.plexos(TORRA2, "TORRA2"),
                 write.to.plexos(TORRA3, "TORRA3"),
                 write.to.plexos(TORRA4, "TORRA4"),
                 write.to.plexos(TORRB1, "TORRB1"),
                 write.to.plexos(TORRB2, "TORRB2"),
                 write.to.plexos(TORRB3, "TORRB3"),
                 write.to.plexos(TORRB4, "TORRB4"),
                 write.to.plexos(NPS, "NPS"),
                 write.to.plexos(WKIEWA1, "WKIEWA1"),
                 write.to.plexos(WKIEWA2, "WKIEWA2"),
                 write.to.plexos(KAREEYA1, "KAREEYA1"),
                 write.to.plexos(KAREEYA2, "KAREEYA2"),
                 write.to.plexos(KAREEYA3, "KAREEYA3"),
                 write.to.plexos(KAREEYA4, "KAREEYA4"),
                 write.to.plexos(HUMENSW, "HUMENSW"),
                 write.to.plexos(HUMEV, "HUMEV"),
                 write.to.plexos(BARRON1, "BARRON-1"),
                 write.to.plexos(BARRON2, "BARRON-2"),
                 write.to.plexos(LK_ECHO, "LK_ECHO"),
                 write.to.plexos(MEADOWBK, "MEADOWBK"),
                 write.to.plexos(FISHER, "FISHER"),
                 write.to.plexos(EILDON1, "EILDON1"),
                 write.to.plexos(EILDON2, "EILDON2"),
                 write.to.plexos(DEVILS_G, "DEVILS_G"),
                 write.to.plexos(BLOWERNG, "BLOWERNG"),
                 write.to.plexos(BASTYAN, "BASTYAN"),
                 write.to.plexos(MACKNTSH, "MACKNTSH"),
                 write.to.plexos(TRIBUTE, "TRIBUTE"),
                 write.to.plexos(CETHANA, "CETHANA"),
                 write.to.plexos(REECE1, "REECE1"),
                 write.to.plexos(REECE2, "REECE2"),
                 write.to.plexos(JBUTTERS, "JBUTTERS"),
                 write.to.plexos(DARTM1, "DARTM1"),
                 write.to.plexos(LD01, "LD01"),
                 write.to.plexos(LD02, "LD02"),
                 write.to.plexos(LD03, "LD03"),
                 write.to.plexos(LD04, "LD04"),
                 write.to.plexos(PTSTAN1, "PTSTAN1"),
                 write.to.plexos(POR03, "POR03"),
                 write.to.plexos(QPS1, "QPS1"),
                 write.to.plexos(QPS2, "QPS2"),
                 write.to.plexos(QPS3, "QPS3"),
                 write.to.plexos(QPS4, "QPS4"),
                 write.to.plexos(MACKAYGT, "MACKAYGT"),
                 write.to.plexos(BBTHREE1, "BBTHREE1"),
                 write.to.plexos(BBTHREE2, "BBTHREE2"),
                 write.to.plexos(BBTHREE3, "BBTHREE3"),
                 write.to.plexos(LADBROK1, "LADBROK1"),
                 write.to.plexos(LADBROK2, "LADBROK2"),
                 write.to.plexos(ROMA_8, "ROMA_8"),
                 write.to.plexos(BDL01, "BDL01"),
                 write.to.plexos(BDL02, "BDL02"),
                 write.to.plexos(VPGS1, "VPGS1"),
                 write.to.plexos(VPGS2, "VPGS2"),
                 write.to.plexos(VPGS3, "VPGS3"),
                 write.to.plexos(VPGS4, "VPGS4"),
                 write.to.plexos(VPGS5, "VPGS5"),
                 write.to.plexos(VPGS6, "VPGS6"),
                 write.to.plexos(JLA01, "JLA01"),
                 write.to.plexos(JLA02, "JLA02"),
                 write.to.plexos(JLA03, "JLA03"),
                 write.to.plexos(JLA04, "JLA04"),
                 write.to.plexos(DRYCGT1, "DRYCGT1"),
                 write.to.plexos(DRYCGT2, "DRYCGT2"),
                 write.to.plexos(DRYCGT3, "DRYCGT3"),
                 write.to.plexos(TVPP104, "TVPP104"),
                 write.to.plexos(JLB01, "JLB01"),
                 write.to.plexos(JLB02, "JLB02"),
                 write.to.plexos(JLB03, "JLB03"),
                 write.to.plexos(YABULU2, "YABULU2"),
                 write.to.plexos(MINTARO, "MINTARO"),
                 write.to.plexos(QPS5, "QPS5"),
                 write.to.plexos(MSTUART3, "MSTUART3"),
                 write.to.plexos(OAKEY1, "OAKEY1"),
                 write.to.plexos(OAKEY2, "OAKEY2"),
                 write.to.plexos(MSTUART1, "MSTUART1"),
                 write.to.plexos(MSTUART2, "MSTUART2"),
                 write.to.plexos(LNGS1, "LNGS1"),
                 write.to.plexos(LNGS2, "LNGS2"),
                 write.to.plexos(YABULU, "YABULU"),
                 write.to.plexos(URANQ11, "URANQ11"),
                 write.to.plexos(URANQ12, "URANQ12"),
                 write.to.plexos(URANQ13, "URANQ13"),
                 write.to.plexos(URANQ14, "URANQ14"),
                 write.to.plexos(BRAEMAR1, "BRAEMAR1"),
                 write.to.plexos(BRAEMAR2, "BRAEMAR2"),
                 write.to.plexos(BRAEMAR3, "BRAEMAR3"),
                 write.to.plexos(BRAEMAR5, "BRAEMAR5"),
                 write.to.plexos(BRAEMAR6, "BRAEMAR6"),
                 write.to.plexos(BRAEMAR7, "BRAEMAR7"),
                 write.to.plexos(CG1, "CG1"),
                 write.to.plexos(CG2, "CG2"),
                 write.to.plexos(CG3, "CG3"),
                 write.to.plexos(CG4, "CG4"),
                 write.to.plexos(MORTLK11, "MORTLK11"),
                 write.to.plexos(MORTLK12, "MORTLK12"),
                 write.to.plexos(TALWA1, "TALWA1"),
                 write.to.plexos(AGLHAL, "AGLHAL"),
                 write.to.plexos(AGLSOM, "AGLSOM"),
                 write.to.plexos(ANGAST1, "ANGAST1"),
                 write.to.plexos(GORDON, "GORDON"),
                 write.to.plexos(GUTHEGA, "GUTHEGA"),
                 write.to.plexos(HVGTS, "HVGTS"),
                 write.to.plexos(LI_WY_CA, "LI_WY_CA"),
                 write.to.plexos(LONSDALE, "LONSDALE"),
                 write.to.plexos(MCKAY1, "MCKAY1"),
                 write.to.plexos(MURRAY, "MURRAY"),
                 write.to.plexos(POAT110, "POAT110"),
                 write.to.plexos(POAT220, "POAT220"),
                 write.to.plexos(POR01, "POR01"),
                 write.to.plexos(SNUG1, "SNUG1"),
                 write.to.plexos(TARRALEA, "TARRALEA"),
                 write.to.plexos(TREVALLN, "TREVALLN"),
                 write.to.plexos(TUMUT3, "TUMUT3"),
                 write.to.plexos(TUNGATIN, "TUNGATIN"),
                 write.to.plexos(UPPTUMUT, "UPPTUMUT"))


fullset[,"Value"] <- as.numeric(as.character(fullset[,"Value"]))