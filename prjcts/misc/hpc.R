################################################
#
#	interface with HPC system
#
################################################

cg.cmd.hpcwrapper<- function(cmd, hpc.walltime=24, hpc.mem='600mb', hpc.nproc=1, hpc.q='pqeph', hpc.load='module load R/2.15')
{
	wrap<- "#!/bin/sh"	
	tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
	wrap<- paste(wrap, tmp, sep='\n')		
	tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
	wrap<- paste(wrap, tmp, sep='\n')
	wrap<- paste(wrap, "#PBS -j oe", sep='\n')
	if(!is.na(hpc.q))
		wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
	wrap<- paste(wrap, hpc.load, sep='\n')
	
	cmd	<- paste(wrap,cmd,sep='\n')
	cmd	
}

cg.cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',outfile,'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',outfile,'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")		
	}
	
}

cg.pipeline.sim.mM<- function()
{
	dummy<- sapply(7000:10000, function(n)
	#dummy<- sapply(10000:19997, function(n)
		{
			cmd		<- paste(HOME,'/startme.R', sep='')
			cmd		<- paste(cmd, ' -exe=MOM.SDE.SIMMO -i=',i,sep='')	
			cmd		<- cg.cmd.hpcwrapper(cmd, hpc.walltime=0, hpc.mem='470mb', hpc.nproc=1, hpc.q='pqeph')	
			file	<- paste('cgsm',paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')	
			cg.cmd.hpccaller( SCRIPT, file, cmd )	
		})				
}

cg.pipeline.pseudodata<- function()
{
	dummy<- sapply(1:12, function(i)
			{
				sapply(1:100, function(n)
						{
							cmd		<- paste(HOME,'/startme.R', sep='')
							cmd		<- paste(cmd, ' -exe=MOM.SDE.PSEUDODATA -i=',i,' -n=',n, sep='')	
							cmd		<- cg.cmd.hpcwrapper(cmd, hpc.walltime=0, hpc.mem='470mb', hpc.nproc=1, hpc.q='pqeph')	
							file	<- paste('cgpd',paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')	
							cg.cmd.hpccaller( SCRIPT, file, cmd )	
						})				
			})	
}

cg.pipeline<- function()
{
	#cg.pipeline.pseudodata()
	cg.pipeline.sim.mM()
}
