library("readxl")
library("tidyverse")

## Load risk assumptions
#BRCA1_Risks <- read_excel('250516-Penetrance_Assumptions.xlsx',1)
#BRCA2_Risks <- read_excel('250516-Penetrance_Assumptions.xlsx',2)
BRCA1_Risks <- read_excel('250520-Penetrance_Assumptions.xlsx',1)
BRCA2_Risks <- read_excel('250520-Penetrance_Assumptions.xlsx',2)


## function to estimated expected risk
fn.expected<-function(cumh_br1, cumh_br2, startage, endage){
    
    mycumh_brca1b<-cumh_br1[endage] -
    cumh_br1[startage+2] +
    (cumh_br1[startage+2] -     cumh_br1[startage+1])/2

    mycumh_brca2b<-cumh_br2[endage] -
    cumh_br2[startage+2] +
    (cumh_br2[startage+2] -     cumh_br2[startage+1])/2

    myout<-cbind(startage, endage, mycumh_brca1b, mycumh_brca2b)
    colnames(myout)<-c("age", "censorage", "brca1H", "brca2H")

    myout
    
}

### Expected 8y risk, chen assumption
myexpected_8<-fn.expected(BRCA1_Risks$"Chen H correction",
            BRCA2_Risks$"Chen H correction",
            30:49,
            pmin(52, (30:49)+8))
## 9y risk etc
myexpected_9<-fn.expected(BRCA1_Risks$"Chen H correction",
            BRCA2_Risks$"Chen H correction",
            30:49,
            pmin(52, (30:49)+9))

myexpected_10<-fn.expected(BRCA1_Risks$"Chen H correction",
            BRCA2_Risks$"Chen H correction",
            30:49,
            pmin(52, (30:49)+10))

myexpected_11<-fn.expected(BRCA1_Risks$"Chen H correction",
            BRCA2_Risks$"Chen H correction",
            30:49,
            pmin(52, (30:49)+11))

myexpected_12<-fn.expected(BRCA1_Risks$"Chen H correction",
            BRCA2_Risks$"Chen H correction",
            30:49,
            pmin(52, (30:49)+12))

myexpected_m<-fn.expected(BRCA1_Risks$"Chen H correction",
            BRCA2_Risks$"Chen H correction",
            30:49,
            rep(52, 20))
##==
### Expected 8y risk, Boadicea assumption
myexpected_8_2<-fn.expected(BRCA1_Risks$"BOADICEA H",
            BRCA2_Risks$"BOADICEA H",
            30:49,
            pmin(52, (30:49)+8))
myexpected_9_2<-fn.expected(BRCA1_Risks$"BOADICEA H",
            BRCA2_Risks$"BOADICEA H",
            30:49,
            pmin(52, (30:49)+9))
myexpected_10_2<-fn.expected(BRCA1_Risks$"BOADICEA H",
            BRCA2_Risks$"BOADICEA H",
            30:49,
            pmin(52, (30:49)+10))
myexpected_11_2<-fn.expected(BRCA1_Risks$"BOADICEA H",
            BRCA2_Risks$"BOADICEA H",
            30:49,
            pmin(52, (30:49)+11))
myexpected_12_2<-fn.expected(BRCA1_Risks$"BOADICEA H",
            BRCA2_Risks$"BOADICEA H",
            30:49,
            pmin(52, (30:49)+12))
myexpected_m_2<-fn.expected(BRCA1_Risks$"BOADICEA H",
            BRCA2_Risks$"BOADICEA H",
            30:49,
            rep(52, 20))
##################

##load protector age/brca distribution file
myagesum<-read_csv(file="250617-protector-agesum.csv")

##################

## Expected events - use assumption 50:50 BRCA1/2 or actual distribution to date from protector
fn.expectevents<-function(insumtab){
    ##doesnt change anything but remove NAs from brca prob
    insumtab[is.na(insumtab)]<-0.5
    
    ## expected events
    ##RRES
    ## assume 50%
    rres1<-insumtab$brca1H %*% insumtab$prres * 0.5 +
        insumtab$brca2H %*% insumtab$prres * 0.5
## use actual distribution
    rres2<-insumtab$brca1H %*% (insumtab$prres * insumtab$prresbrca1)+
        insumtab$brca2H %*% (insumtab$prres * (1-insumtab$prresbrca1))
##RRSO
## assume 50%
    rrso1<-insumtab$brca1H %*% insumtab$prrso * 0.5 +
        insumtab$brca2H %*% insumtab$prrso * 0.5
## use actual distribution
    rrso2<-insumtab$brca1H %*% (insumtab$prrso * insumtab$prrsobrca1)+
        insumtab$brca2H %*% (insumtab$prrso * (1-insumtab$prrsobrca1))
    c(rres1, rres2, rrso1, rrso2)
    }

## summary of assumptions
myscenfu<-rbind(
    fn.expectevents(cbind(myagesum, myexpected_8[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_9[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_10[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_11[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_12[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_m[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_8_2[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_9_2[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_10_2[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_11_2[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_12_2[,c(2,3,4)])),
    fn.expectevents(cbind(myagesum, myexpected_m_2[,c(2,3,4)]))
    )

## write file
write.csv(myscenfu[,c(2,4)], file="250516-riskscenarios.csv")

## generate summary table for supplement
mysumtab8<- cbind(myagesum, myexpected_8[,c(2,3,4)])
mysumtab8_2<- cbind(myagesum, myexpected_8_2[,c(2,3,4)])

mysumtab10<- cbind(myagesum, myexpected_10[,c(2,3,4)])
mysumtab10_2<- cbind(myagesum, myexpected_10_2[,c(2,3,4)])

mysumtab_all<-cbind(mysumtab8,
                    mysumtab8_2[,7:8],
                    mysumtab10[,7:8],
                    mysumtab10_2[,7:8])
colnames(mysumtab_all)[7:14]<-c("brca1H_8_A", "brca2H_8_A", "brca1H_8_B", "brca2H_8_B", "brca1H_10_A", "brca2H_10_A", "brca1H_10_B", "brca2H_10_B")

mysumtab_all<-mysumtab_all[,c(1:6,7,11,8,12,9,13,10,14)]

    fn.format<-function(ind, ndigit=2){                                                       format(round(ind,ndigit), nsmall=ndigit)
    }

mysumtab_all[,c(2,3,4,5,7:14)]<-
    fn.format(100*mysumtab_all[,c(2,3,4,5,7:14)],2)

write.csv(mysumtab_all, file="250517-riskassumptions.csv")





