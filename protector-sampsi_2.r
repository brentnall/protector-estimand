## Sample size calcs for scenarios

## Libraries
library(future.apply)
library("binom")
library("tidyverse")
## Set up workers to run calculations in parallel
plan(multisession, workers = 7)  

## functions

##format dps
fn.format<-function(ind, ndigit=2){                                                       format(round(ind,ndigit), nsmall=ndigit)
}

## Ratio cancers prevented
mytstat<-function(myorres, myorrso, myerres, myerrso){
    (1-(myorres/myerres)) /     (1-(myorrso / myerrso))
    }

## bootstrap CI for ratio cancers prevented
fn.bootstat<-function(inn, inp, ine, inR=1000, inci=c(0.025, 0.975)){

    myo1<-rbinom(inR, inn[1], inp[1])
    myo2<-rbinom(inR, inn[2], inp[2])

    mybtstat <- mytstat(myo1, myo2, rep(ine[1],inR), rep(ine[2], inR))

    myci<-quantile(mybtstat,inci)  

    }


## power for o/e main anal + thresholds 80%/90% power non-inferiority
fn.powercalc1<-function(inscenario){

    ##debug
##    inscenario<-myscenario
    
    myn<-inscenario$ntot * (1-inscenario$dropout)

    mye1<-inscenario$e1
    mye2<-inscenario$e2

    myeff1<-inscenario$eff1
    myeff2<-inscenario$eff2

    myprres<-inscenario$prres
    myprrso<-inscenario$prrso
    myprest<-1-myprres-myprrso

    mynrres<-round(myn*myprres,0)
    mynrrso<-round(myn*myprrso,0)
    
    myB<-inscenario$B

    myorres<-rbinom(myB, mynrres, mye1*(1-myeff1))
    myorrso<-rbinom(myB, mynrrso, mye2*(1-myeff2))

    myerres<-mye1*mynrres
    myerrso<-mye2*mynrrso

    myci<-sapply(0:mynrres, function(idx) binom.wilson(idx,mynrres)$upper)
    mycit<-myci*mynrres/myerres <0.8
    mydpois<-dpois(0:mynrres, myerres*(1-myeff1))
    outpower1<-mydpois%*%mycit

    myrat<-mytstat(myorres, myorrso, myerres, myerrso)

    mybci <- future_sapply(1:myB,
                       function(idx) fn.bootstat(
                                         c(mynrres,mynrrso),
                                         c(myorres[idx]/mynrres, myorrso[idx]/mynrrso),
                                         c(myerres, myerrso),
                                         inscenario$R),
                       future.seed=TRUE)

############## thresholds
    outthresh2<-quantile(mybci[1,], c(0.1,0.2)) ## 90% power, 80% power
    list(outpower1, outthresh2)

    }

#########
## Scenarios
## original
myscenario<-list(e1=0.029,
                 e2=0.044,
                 eff1=0.65,
                 eff2=0.96,
                 prres=0.454,
                 prrso=0.286,
                 dropout=0.05,
                 ntot=2533,
                 B=10000,
                 R=1000)

## original calc from first draft
myscen1<-fn.powercalc1(myscenario)
myscen1

##################
## systematic exploration of scenarios
myscenarios<-read.csv("250516-riskscenarios.csv")

## tables varying risk assumptions
fn.varyrisk<-function(myscenario){
    
    myscen_all<-list()
    myscenario$eff1<-0.65
    for(idx in 1:nrow(myscenarios)){
        myscenario$e1<-myscenarios[idx,2]
        myscenario$e2<-myscenarios[idx,3]
        myscen_all[[idx]]<-fn.powercalc1(myscenario)
    }
    
    ## 60% efficacy
    myscen_all2<-list()
    myscenario$eff1<-0.6
    for(idx in 1:nrow(myscenarios)){
        myscenario$e1<-myscenarios[idx,2]
        myscenario$e2<-myscenarios[idx,3]
        myscen_all2[[idx]]<-fn.powercalc1(myscenario)
    }
    
    ##55% efficacy
    myscen_all3<-list()
    myscenario$eff1<-0.55
    for(idx in 1:nrow(myscenarios)){
        myscenario$e1<-myscenarios[idx,2]
        myscenario$e2<-myscenarios[idx,3]
        myscen_all3[[idx]]<-fn.powercalc1(myscenario)
    }

   
    mypower_all<-cbind(unlist(myscen_all),
                       unlist(myscen_all2),
                       unlist(myscen_all3))
    

    mypower<-mypower_all[seq(1,nrow(mypower_all),by=3),]
    powersens<-data.frame(cbind(1:12, fn.format(100*myscenarios$V1,1), fn.format(100*mypower,0)))
    colnames(powersens)<-c("Scenario", "Risk", "Power (A)", "Power (B)", "Power (C)")

    ## secondary info, non-inf threshold 90% power
    mynoninf<-mypower_all[seq(2,nrow(mypower_all),by=3),]
    noninf<-data.frame(cbind(1:12, fn.format(100*myscenarios$V1,1), fn.format(myscenarios$V2*100,1), fn.format(mynoninf*100,0)))
    colnames(noninf)<-c("Scenario", "Risk RRES", "Risk RRSO", "Threshold (A)", "Threshold (B)", "Threshold (C)")

    list(powersens, noninf)
}

## run scenarios
## n=1000 for rres
myscenario$ntot=round(1000/0.454)+1
myoutA<-fn.varyrisk(myscenario)
## n=1150 for rres
myscenario$ntot=round(1150/0.454)+1
myoutB<-fn.varyrisk(myscenario)
## n=1300 for rres
myscenario$ntot=round(1300/0.454)+1
myoutC<-fn.varyrisk(myscenario)

## summary table
mysumsamp<-data.frame(
    cbind(
        rep(c(1000,1150,1300), each=12),
        rbind(
            myoutA[[1]],
            myoutB[[1]],
    myoutC[[1]]
    )))

colnames(mysumsamp)[1]<-c("Sample size")

write_csv(mysumsamp %>% arrange(Risk,`Sample size`, Scenario ), "250515-risktable-sampsize-power.csv")
    
## summary table 2
mysumsamp2<-data.frame(
    cbind(
        rep(c(1000,1150,1300), each=12),
        rbind(
            myoutA[[2]],
            myoutB[[2]],
    myoutC[[2]]
    )))

colnames(mysumsamp2)[1]<-c("Sample size")

write_csv(mysumsamp2 %>% arrange(Risk.RRES,`Sample size`, Scenario ), "250515-risktable-sampsize-threshold.csv")
    
