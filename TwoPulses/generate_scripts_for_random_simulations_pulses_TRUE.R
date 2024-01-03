##This is the main script for running the simulations for training the neural network for the Two Pulses model
##Software and packages needed:
## (i) Slim as a module of the cluster
## (ii) R and packages ggplot2, gridExtra and truncnorm 
#### On HPC runs on anaconda3/personal. Need to run: module load anaconda3/personal; anaconda- setup; conda install r; conda install r-ggplot2; conda install r-gridextra.

library(ggplot2)
library(gridExtra)
library(truncnorm)

##Define Path
path="./"

slim_model="SLIM_AutX_pulses.slim"
totalpop <- c(1000)
CLM=c(0.082,0.277)
MXL=c(0.043,0.500)
PEL=c(0.032,0.761)
PUR=c(0.138, 0.150) #AFR,NAT
ACB=c(0.881,0.005)
ASW=c(0.765,0.041)
pops=cbind(CLM,MXL,PEL,PUR,ACB,ASW)

GEN <- 19 #generations
ITE <- 10000 # runs with different seeds
LA0 <- 2873038360
LX0 <- 3028738110
scale <- c("1000")
LAstep <- 100


lastrun_pop=c(0,0,0,0)
for (pop in 1:length(colnames(pops))){
	if (!(dir.exists(paste(path,colnames(pops)[pop],sep="")))){
		dir.create(paste(path,colnames(pops)[pop],sep=""))
		dir.create(paste(path,colnames(pops)[pop],"/ancestries",sep=""))
		dir.create(paste(path,colnames(pops)[pop],"/fragments",sep=""))		
	}
	if (!(dir.exists(paste(pathEph,colnames(pops)[pop],sep="")))){
		dir.create(paste(pathEph,colnames(pops)[pop],sep=""))
		#dir.create(paste(pathEph,colnames(pops)[pop],"/ancestries",sep=""))
		#dir.create(paste(pathEph,colnames(pops)[pop],"/fragments",sep=""))
	}
  AMcomb_pop_file<-paste(path,"AM_combinations_pulses_true_scale_",scale,"_pop_",pop,".txt",sep="")
	if (file.exists(AMcomb_pop_file)){
		AMcomb_pop=read.table(AMcomb_pop_file)
		lastrun_pop[pop]=max(AMcomb_pop$V1)
	}
}

system(paste("rm ",path,"runrunfile_slimAutX_pulses_true_scale_10E3.sh",sep=""))
for (SCALE in scale) {
for (TOTALPOP in totalpop) {
for (pop in 1:length(colnames(pops))){  
      RATPOP1=pops[1,pop] #NAT in PUR
      RATPOP2=pops[2,pop] #AFR in PUR
      POP1 <- round(TOTALPOP*RATPOP1)
      POP2 <- round(TOTALPOP*RATPOP2)# population
      POP3 <- round(TOTALPOP-POP1-POP2)
      #Name input and output files, define arguments for slim and set up the running script
      output_jobfile <- paste(path,colnames(pops)[pop],"/Slim_output_AutX_pulses_true_","_GEN_",as.character(GEN),"_POP1_",as.character(RATPOP1),"_POP2_",as.character(RATPOP2),"_totalpop_",as.character(TOTALPOP),"_pop_",pop,"_",colnames(pops)[pop],"_scale_",SCALE,"_run_$(($PBS_ARRAY_INDEX+",lastrun_pop[pop],")).txt",sep="")
      fragsfile <- paste(path,colnames(pops)[pop],"/fragments/Frags_AutX_pulses_true_","_GEN_",as.character(GEN),"_POP1_",as.character(RATPOP1),"_POP2_",as.character(RATPOP2),"_totalpop_",as.character(TOTALPOP),"_pop_",pop,"_",colnames(pops)[pop],"_scale_",SCALE,"_run_$(($PBS_ARRAY_INDEX+",lastrun_pop[pop],")).txt",sep="")
      ancfile <- paste(path,colnames(pops)[pop],"/ancestries/Ancestries_AutX_pulses_true_","_GEN_",as.character(GEN),"_POP1_",as.character(RATPOP1),"_POP2_",as.character(RATPOP2),"_totalpop_",as.character(TOTALPOP),"_pop_",pop,"_",colnames(pops)[pop],"_scale_",SCALE,"_run_$(($PBS_ARRAY_INDEX+",lastrun_pop[pop],")).txt",sep="")
      
      command_runjobfile <- paste("slim -d AM1_t1=1000",
                                        " -d AM2_t1=1000",
                                        " -d AM3_t1=1000",
                                        " -d SB1_t1=0.0",
                                        " -d SB2_t1=0.0",
				  " -d POP=",pop,
	                          " -d GEN=",GEN,
	                          " -d RUN=$(($PBS_ARRAY_INDEX+",lastrun_pop[pop],"))",
	                          " -d LA0=",LA0,
	                          " -d LX0=",LX0,
	                          " -d SCALE=",SCALE,
	                          " -d LAstep=",LAstep,
	                          " -d totalpop=",TOTALPOP,
	                          " -d pop1_ratio=",RATPOP1,
	                          " -d pop2_ratio=",RATPOP2,
                                        " -d prop_pop1_t3=0.0",
                                        " -d prop_pop2_t3=0.0",
                                        " -d prop_pop3_t3=0.0",
 	                          " ", 		                          	                          
	                          pathModel,slim_model," > ",output_jobfile,sep="")
      command_split_output_1 <- paste("linefrag=$(grep -n ^###Frag ",output_jobfile," | sed -r 's/([0-9]+):.*/\\","1","/g')",sep="")
      command_split_output_2 <- paste("tail -n+$(($linefrag+1)) ",output_jobfile," > ",fragsfile)
      command_split_output_3 <- paste("lineanc=$(grep -n ^###Ancestry ",output_jobfile," | sed -r 's/([0-9]+):.*/\\","1","/g')",sep="")
      command_split_output_4 <- paste("tail -n+$(($lineanc+1)) ",output_jobfile," | head -n $(($linefrag-$lineanc-1)) > ",ancfile)

      
      runjobfile <-paste(path,colnames(pops)[pop],"/run_job_AutX_pulses_true_GEN_",as.character(GEN),"_POP1_",as.character(RATPOP1),"_POP2_",as.character(RATPOP2),"_totalpop_",as.character(TOTALPOP),"_pop_",pop,"_",colnames(pops)[pop],"_scale_",SCALE,".sh",sep="")
      jobname <- paste("job_pulses_true_GEN_",as.character(GEN),"_POP1_",as.character(RATPOP1),"_POP2_",as.character(RATPOP2),"_totalpop_",as.character(TOTALPOP),"_pop_",pop,"_",colnames(pops)[pop],"_scale_",SCALE,sep="")
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
      #wait if 50 jobs
      systemcommandNJOBS <- "NJOBS=$(qstat | grep job | wc -l); while [ $NJOBS -gt 47 ]; do echo sleep; sleep 30; NJOBS=$(qstat | grep job | wc -l); done"
      write(systemcommandNJOBS,paste(path,"runrunfile_slimAutX_pulses_true_scale_10E3.sh",sep=""),append=TRUE)
      write(qsubcommand,paste(path,"runrunfile_slimAutX_pulses_true_scale_10E3.sh",sep=""),append=TRUE)
      system(paste("chmod 770 ",path,"runrunfile_slimAutX_pulses_true_scale_10E3.sh",sep=""))
}}}



