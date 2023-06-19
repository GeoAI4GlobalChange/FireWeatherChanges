rm(list=ls())
run_cases=c(1,2)
target_region = 'all' 
for (case in run_cases)
{
rm(list=ls()[! ls() %in% c("case","run_cases","target_region")])#rm(list=ls())

beta_values=array(0,c(1,3))

dir='direcory_to_the_code/'
source(paste(dir,'ECOF.r',sep=""))
ls()
dir_data='direcory_to_the_data/'
dir_beta='direcory_to_be_saved/'
if (case==1){
Z<-readin(paste(dir_data,'sig_ltn_obs_all_',target_region,'.txt',sep=""),paste(dir_data,'noise1_ltn_',target_region,'.txt',sep=""),paste(dir_data,'noise2_ltn_',target_region,'.txt',sep=""))
}
else if (case==2){
Z<-readin(paste(dir_data,'sig_human_obs_all_',target_region,'.txt',sep=""),paste(dir_data,'noise1_human_',target_region,'.txt',sep=""),paste(dir_data,'noise2_human_',target_region,'.txt',sep=""))
}


############################################################

Zr=Z
nsig =1
o1.rof<-tls.ROF(Zr@Y,Zr@X,Zr@noise1,Zr@noise2,nsig=nsig,plev=0.90,nsim.CI=1000,nsim.rcc=1000)

print(o1.rof$beta.CI)

beta_values[,]=o1.rof$beta.CI


if (case==1){
saveRDS(beta_values,file=paste(dir_beta,'beta_values_ltn_all','_',target_region,'.rds',sep=''))
}
else {
saveRDS(beta_values,file=paste(dir_beta,'beta_values_human_all','_',target_region,'.rds',sep=''))
}

}
