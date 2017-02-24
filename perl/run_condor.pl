#!/usr/bin/perl -w
#######################################################
# run R script on the cluster                         #
#######################################################
my $HOME="/users/ecologie/camacho";
my $projectDir="$HOME/work/projects/ReinfectionProbability";
my $Condor_job_identifier="Rscript";

#------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------------#

my @model_list= ("AoN","PPI","SIRS");
# my @model_list= ("SIRS");
my $analysis= "large";

#------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------------#

my $scriptDir="$projectDir/R";
my $Rscript= "dev.r";

# my @node_list= ("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20");
my @node_list= ("19");

my $exe= "$HOME/bin/R";
my $uni= "vanilla";
my $mem= "1000";
my $replicate= "4";

#request reserved machines
my $my_requirements="";
foreach $node (@node_list)
{
  $my_requirements .= "(machine == \"ecoevo$node.ecoevo.biologie.ens.fr\")+";
}
#remove the + at the end of the string
$my_requirements=substr($my_requirements,0,-1);

foreach $model (@model_list)
{

  my $jobDir="$projectDir/cluster_output/$model";

  # remove the .txt (submission, out, log and error files)
  system("mkdir -p $jobDir");
  system("rm -r $jobDir/*.txt");


  #my $scriptDir= "$jobDir";
  ## system("rm -r $jobDir/Rsave $jobDir/pdf"); 
  # system("mkdir -p $jobDir");
  my $in= "$scriptDir/$Rscript";
  my $out= "$jobDir/out_R\$(Process).txt";

  my $Condor_jobname="$Condor_job_identifier\_$model";
  print "generate $Condor_jobname\n";
  open(Condor_SCRIPT,">$Condor_jobname.txt");
  print Condor_SCRIPT "executable = $exe\n";
  print Condor_SCRIPT "universe = $uni\n";
  print Condor_SCRIPT "request_memory = $mem\n";
  print Condor_SCRIPT "+RequiresWholeMachine = False\n";
    # print Condor_SCRIPT "+RequiresWholeMachine = True\n";
    print Condor_SCRIPT "Requirements = $my_requirements\n";
    #print Condor_SCRIPT "Requirements = machine == \"ecoevo13.ecoevo.biologie.ens.fr\"\n";
    print Condor_SCRIPT "getenv = True\n";

    print Condor_SCRIPT "arguments = --vanilla\n";
	#print Condor_SCRIPT "transfer_files = ALWAYS\n";
    #print Condor_SCRIPT "transfer_input_files = $in\n";
    print Condor_SCRIPT "notification = NEVER\n";
    print Condor_SCRIPT "input = $in\n";
    print Condor_SCRIPT "output = $out\n";
    print Condor_SCRIPT "error = $jobDir/error_R\$(Process).txt\n";
    print Condor_SCRIPT "log = $jobDir/log_R\$(Process).txt\n";
    print Condor_SCRIPT "environment = process=\$(Process);replicate=$replicate;model=$model;analysis=$analysis;jobDir=$jobDir\n";
    # print Condor_SCRIPT "environment = analysis=$analysis\n";


    print Condor_SCRIPT "queue $replicate\n";
    
    print "submit $Condor_jobname with $replicate replicate\n";
    
    system("condor_submit $Condor_jobname.txt")==0
    or die "failed submitting Condor job \n";
    
    #### GIVE THE JOB A SECOND TO GET SUBMITTED ####
    sleep(1);
    print "submitted.\n";
  }

