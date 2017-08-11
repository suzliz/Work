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
        P_MW <- Av_MW * (1-PDR)
        
        ## Calculate probability of no outage
        prob_avail <- (1-prob_out)*(1-prob_par)
        ## Calculate probability of a partial outage
        prob_pout <- (1-prob_out)*prob_par
        
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
        Comb_Prob <- round(Comb_Prob,5)
        
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
        names(Comb) <- c("Available", "Partial", "Full_Out","Total_Cap","Out_Fac","MTTR","FOR")
        
        ## Only return probabilities above zero (According to specified significance)
        if (rep_zero==FALSE)
                return(Comb[Comb[,"FOR"]> 0,])
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
fullexp_prob <- function(unit_cap,n,prob_out,prob_par,PDR,MTTR_O,MTTR_P,MTTR_D=0.5,CCGT=FALSE){
        prob_avail <- (1-prob_out)*(1-prob_par)
        prob_pout <- (1-prob_out)*prob_par
        prob <- c(prob_avail, prob_pout,prob_out)
        
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
        Comb_Prob <- round(Comb_Prob,5)
        
        ## Combine all data
        
        ## Calculate Percentage of Capacity on Outage
        CapOut <- round((sum(unit_cap) - Comb_Cap)/sum(unit_cap)*100,2)
        Comb <- data.frame(cbind(Comb_State,Comb_Cap,CapOut, Comb_Prob,MTTR))
        ## Sum up probabilities of equal operating states to reduce overall size
        ## i.e. where the same configuration of operating states has the same available capacity
        Reduced <- aggregate(Comb$Comb_Prob,by=list(Comb$Avail,Comb$Par,Comb$Out,Comb$Comb_Cap,Comb$CapOut, Comb$MTTR),FUN=sum)
        names(Reduced) <- c("Available", "Partial", "Full_Out","Total_Cap","Out_Fac","MTTR","FOR")
        ## Only return probabilities above zero (5dp)
        Reduced[Reduced[,"FOR"]> 0,]
       
}








Murray <- Outage_Calc(1500,14,0.0082217,0.0000333,0.119697,20.45*2,5.9*2)
Hallet <- Outage_Calc(220,12,0.0066015,0.0008687,0.3189,33.45*2,10.75*2)
SHGEN <- Outage_Calc(240,4,0.0082217,0.0000333,0.119697,20.45*2,5.9*2)

