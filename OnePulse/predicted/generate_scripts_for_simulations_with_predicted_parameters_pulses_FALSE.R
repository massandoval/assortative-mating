##This is the main script for running the simulations for training the neural network for the OnePulse model
##Software and packages needed:
## (i) Slim as a module of the cluster
## (ii) R and packages ggplot2, gridExtra and truncnorm 
#### On HPC runs on anaconda3/personal. Need to run: module load anaconda3/personal; anaconda- setup; conda install r; conda install r-ggplot2; conda install r-gridextra.

library(ggplot2)
library(gridExtra)
library(truncnorm)

##Define Path
path="./"
##Define Base model
slim_model="SLIM_AutX_pulses_predicted.slim"
totalpop <- c(1000)
CLM=c(0.082,0.277)#c(0.089,0.278)
MXL=c(0.043,0.500)#c(0.048,0.502)
PEL=c(0.032,0.761)#c(0.033,0.768)
PUR=c(0.138, 0.150)#c(0.152,0.147) #AFR,NAT
ACB=c(0.881,0.005)
ASW=c(0.765,0.041)
pops=cbind(CLM,MXL,PEL,PUR,ACB,ASW)
POP_names=c("CLM","MXL","PEL","PUR","ACB","ASW")
GEN <- 19 #generations
ITE <- 25 # runs with different seeds
LA0 <- 2873038360
LX0 <- 3028738110
SCALE <- "1000"
TOTALPOP<-1000
LAstep <- 100
system(paste("rm ",path,"runrunfile_slimAutX_pulses_false_predicted_scale_10E3.sh",sep=""))
for (pop in 1:length(colnames(pops))){
  RATPOP1=pops[1,pop] #NAT in PUR
  RATPOP2=pops[2,pop] #AFR in PUR
  POP1 <- round(TOTALPOP*RATPOP1)
  POP2 <- round(TOTALPOP*RATPOP2)# population
  POP3 <- round(TOTALPOP-POP1-POP2)
  param_predicted<-read.table("./mean_parameter_allpop_false.txt",header=TRUE)
	POPname=POP_names[pop]
	AM1=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="AM1")&(param_predicted[,"pop"]==POPname)],digits=5))
	AM2=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="AM2")&(param_predicted[,"pop"]==POPname)],digits=5))
	AM3=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="AM3")&(param_predicted[,"pop"]==POPname)],digits=5))
	SB1=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="SB1")&(param_predicted[,"pop"]==POPname)],digits=5))
	SB2=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="SB2")&(param_predicted[,"pop"]==POPname)],digits=5))
	prop_pop1_t1=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="prop_pop1_t1")&(param_predicted[,"pop"]==POPname)],digits=5))
	prop_pop2_t1=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="prop_pop2_t1")&(param_predicted[,"pop"]==POPname)],digits=5))
	prop_pop3_t1=as.character(round(param_predicted[,"value"][(param_predicted[,"parameter"]=="prop_pop3_t1")&(param_predicted[,"pop"]==POPname)],digits=5))
	print(paste("pop",POPname,sep=":"))
	print(paste("AM1",AM1,sep=":"))
	print(paste("AM2",AM2,sep=":"))
	print(paste("AM3",AM3,sep=":"))
	print(paste("SB1",SB1,sep=":"))
 	print(paste("SB2",SB2,sep=":"))
	print(paste("POP1",POP1,sep=":"))
	print(paste("POP2",POP2,sep=":"))
	print(paste("POP3",POP3,sep=":"))    
	output_jobfile <- paste(path,colnames(pops)[pop],"/Slim_output_AutX_pulses_false_predicted_",colnames(pops)[pop],"_AM1_",as.character(AM1),"_AM2_",as.character(AM2),"_AM3_",as.character(AM3),"_SB1_",as.character(SB1),"_SB2_",as.character(SB2),"_GEN_",as.character(GEN),"_POP1_",as.character(POP1),"_POP2_",as.character(POP2),"_POP3_",as.character(POP3),"_scale_",SCALE,"_run_${PBS_ARRAY_INDEX}.txt",sep="")
	fragsfile <- paste(path,colnames(pops)[pop],"/fragments/Frags_AutX_pulses_false_predicted_",colnames(pops)[pop],"_AM1_",as.character(AM1),"_AM2_",as.character(AM2),"_AM3_",as.character(AM3),"_SB1_",as.character(SB1),"_SB2_",as.character(SB2),"_GEN_",as.character(GEN),"_POP1_",as.character(POP1),"_POP2_",as.character(POP2),"_POP3_",as.character(POP3),"_scale_",SCALE,"_run_${PBS_ARRAY_INDEX}.txt",sep="")
	ancfile <- paste(path,colnames(pops)[pop],"/ancestries/Ancestries_AutX_pulses_false_predicted_",colnames(pops)[pop],"_AM1_",as.character(AM1),"_AM2_",as.character(AM2),"_AM3_",as.character(AM3),"_SB1_",as.character(SB1),"_SB2_",as.character(SB2),"_GEN_",as.character(GEN),"_POP1_",as.character(POP1),"_POP2_",as.character(POP2),"_POP3_",as.character(POP3),"_scale_",SCALE,"_run_${PBS_ARRAY_INDEX}.txt",sep="")
	command_runjobfile <- paste("slim -d POP=",pop,
    " -d RUN=$PBS_ARRAY_INDEX",
    " -d AM1_t1=",AM1,
    " -d AM2_t1=",AM2,
    " -d AM3_t1=",AM3,
    " -d SB1_t1=",SB1,
    " -d SB2_t1=",SB2,
    " -d GEN=",GEN,
    " -d ITE=",ITE,
    " -d LA0=",LA0,
    " -d LX0=",LX0,
    " -d SCALE=",SCALE,
    " -d LAstep=",LAstep,
    " -d totalpop=",TOTALPOP,
    " -d pop1_ratio=",RATPOP1,
    " -d pop2_ratio=",RATPOP2,
              " -d prop_pop1_t1=1.0",
              " -d prop_pop2_t1=1.0",
              " -d prop_pop3_t1=1.0",
              " -d prop_pop1_t2=0",
              " -d prop_pop2_t2=0",
              " -d prop_pop3_t2=0",     
              " -d prop_pop1_t3=0",
              " -d prop_pop2_t3=0",
              " -d prop_pop3_t3=0",     
    " ",
    "../../slim_script/",slim_model," > ",output_jobfile,sep="")
	command_split_output_1 <- paste("linefrag=$(grep -n ^###Frag ",output_jobfile," | sed -r 's/([0-9]+):.*/\\","1","/g')",sep="")
	command_split_output_2 <- paste("tail -n+$(($linefrag+1)) ",output_jobfile," > ",fragsfile)
	command_split_output_3 <- paste("lineanc=$(grep -n ^###Ancestry ",output_jobfile," | sed -r 's/([0-9]+):.*/\\","1","/g')",sep="")
	command_split_output_4 <- paste("tail -n+$(($lineanc+1)) ",output_jobfile," | head -n $(($linefrag-$lineanc-1)) > ",ancfile)

	      
	runjobfile <-paste(path,colnames(pops)[pop],"/run_job_pulses_false_AM1_",as.character(AM1),"_AM2_",as.character(AM2),"_AM3_",as.character(AM3),"_SB1_",as.character(SB1),"_SB2_",as.character(SB2),"_GEN_",as.character(GEN),"_POP1_",as.character(POP1),"_POP2_",as.character(POP2),"_POP3_",as.character(POP3),"_scale_",SCALE,".sh",sep="")
	jobname <- paste("job_pulses_false_predicted_",colnames(pops)[pop],"_AM1_",as.character(AM1),"_AM2_",as.character(AM2),"_AM3_",as.character(AM3),"_SB1_",as.character(SB1),"_SB2_",as.character(SB2),"_GEN_",as.character(GEN),"_POP1_",as.character(POP1),"_POP2_",as.character(POP2),"_POP3_",as.character(POP3),sep="")
	write("#PBS -l walltime=71:59:00",runjobfile)
	write("#PBS -l select=1:ncpus=1:mem=10gb",runjobfile,append=TRUE)
  write(paste("#PBS -J 1-",ITE,sep=""),runjobfile,append=TRUE)
	write(paste("cd ",path,sep=""),runjobfile,append=TRUE)
	write("module load slim",runjobfile,append=TRUE)
	write(command_runjobfile,runjobfile,append=TRUE)
	write(command_split_output_1,runjobfile,append=TRUE)
	write(command_split_output_2,runjobfile,append=TRUE)
	write(command_split_output_3,runjobfile,append=TRUE)
	write(command_split_output_4,runjobfile,append=TRUE)      
 	system(paste("chmod 770 ",runjobfile,sep=""))
 	qsubcommand <- paste("qsub -N \"",jobname,"\" ",runjobfile,sep="")
 	systemcommandNJOBS <- "NJOBS=$(qstat | grep job | wc -l); while [ $NJOBS -gt 47 ]; do echo sleep; sleep 30; NJOBS=$(qstat | grep job | wc -l); done"
	write(systemcommandNJOBS,paste(path,"runrunfile_slimAutX_pulses_false_predicted_scale_10E3.sh",sep=""),append=TRUE)
	write(qsubcommand,paste(path,"runrunfile_slimAutX_pulses_false_predicted_scale_10E3.sh",sep=""),append=TRUE)

}



