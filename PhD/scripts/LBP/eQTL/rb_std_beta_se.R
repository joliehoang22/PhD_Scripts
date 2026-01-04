# step1: select the same number of eGenes from both datasets
data1=fread(paste0(INDIRT,data1_name,"_and_",data2_name,".txt"),head=T,stringsAsFactors=F,data.table=F)
index=which(data1$SE!=0 & data1$se2!=0) #return row indices where condition is true 
#SE = standard error column for dataset1 and se2 = standard error column for dataset2 
data1=data1[index,]
index1=which(data1$p<5e-08);
index2=which(data1$p2<5e-08);

shared_index=intersect(index1,index2)
unique_index1=index1[!index1 %in% shared_index]
unique_index2=index2[!index2 %in% shared_index]

if(length(unique_index2)<length(unique_index1)){
    tmp_index=c(shared_index,unique_index2,sample(unique_index1,length(unique_index2)))
}else{
    tmp_index=c(shared_index,unique_index1,sample(unique_index2,length(unique_index1)))
}
data2=data1[tmp_index,] 


# step2: standarized beta and SE in SD unit
# z represent z statistics, p presents allele frequency, n represents sample size.\
calcu_std_b_se<-function(z,p,n){
    std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
    std_se=1/sqrt(2*p*(1-p)*(n+z^2))
    res<-data.frame(std_b_hat,std_se);
    return(res)
}
#this is fine to skipped as long as SEs and betas are not mismatched in scale right? Standardization changes the units, not the correlation.

###################### 
###################### Alternative approach: Select from one reference dataset only
reference_index <- which(data1$p < 5e-08)  # or data1$p2 < 5e-08
data2 <- data1[reference_index, ]
# Then proceed with step 3
###################### 
###################### 

# step3: calcualte rb
# b1 and se1 represent the estimate and SE for eQTLs across probes in one tissue, b2 and se2 represent the estimate and SE for eQTLs across probes in the other tissue
# theta = sample overlap * phnotypic correlation; theta could also be estimated from null SNPs; if two samples were independent, theta = 0.
# Please note that the effect allele of one SNP between two tissues should be the same.
calcu_cor_true<-function(b1,se1,b2,se2,theta){
    idx=which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2));
    if(length(idx)>0){
        b1=b1[-idx];se1=se1[-idx]
        b2=b2[-idx];se2=se2[-idx]
        theta=theta[-idx]
    }

    var_b1=var(b1,na.rm=T)-mean(se1^2,na.rm=T)
    var_b2=var(b2,na.rm=T)-mean(se2^2,na.rm=T)
    if(var_b1<0){
      var_b1=var(b1,na.rm=T)
    }
    if(var_b2<0){
      var_b2=var(b2,na.rm=T)
    }
    cov_b1_b2=cov(b1,b2,use="complete.obs")-mean(theta,na.rm=T)*sqrt(mean(se1^2,na.rm=T)*mean(se2^2,na.rm=T))
    r=cov_b1_b2/sqrt(var_b1*var_b2)

    r_jack=c()
    n=length(b1)
    for(k in 1:n) {
       b1_jack=b1[-k];se1_jack=se1[-k];var_b1_jack=var(b1_jack,na.rm=T)-mean(se1_jack^2,na.rm=T)
       b2_jack=b2[-k];se2_jack=se2[-k];var_b2_jack=var(b2_jack,na.rm=T)-mean(se2_jack^2,na.rm=T)
       if(var_b1_jack<0){
        var_b1_jack=var(b1_jack,na.rm=T);
       }
       if(var_b2_jack<0){
        var_b2_jack=var(b2_jack,na.rm=T);
       }
       theta_jack=theta[-k];
       cov_e1_jack_e2_jack=mean(theta_jack,na.rm=T)*sqrt(mean(se1_jack^2,na.rm=T)*mean(se2_jack^2,na.rm=T))
       cov_b1_b2_jack=cov(b1_jack,b2_jack,use="complete.obs")-cov_e1_jack_e2_jack
       r_tmp=cov_b1_b2_jack/sqrt(var_b1_jack*var_b2_jack)
       r_jack=c(r_jack,r_tmp)
    }
    r_mean=mean(r_jack,na.rm=T)
    idx=which(is.na(r_jack))
    if(length(idx)>0){
        se_r=sqrt((n-1)/n*sum((r_jack[-idx]-r_mean)^2))
    }else{
    	se_r=sqrt((n-1)/n*sum((r_jack-r_mean)^2))
    }
    res<-cbind(r,se_r)
    return(res)
}



