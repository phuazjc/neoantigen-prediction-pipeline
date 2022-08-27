wd <- "/mnt/projects/zhaiww1/planet/tmp4cheryl/test/"

bam.folder <- "/mnt/projects/zhaiww1/planet/tmp4cheryl/Liver_DNA_bam/"
mutect.folder <- "/mnt/projects/zhaiww1/planet/tmp4cheryl/Liver_DNA_bam/mutect"
strelka.folder <- "/mnt/projects/zhaiww1/planet/DNA/liver_tcr_hch/all_sequenced_so_far_20180412/strelka_results"

# Selected samples
samples <- c("WHT165","WHT166","WHT167","WHT169","WHT170")

# Read in bam files
bamfiles = list.files(path=bam.folder, pattern="^normal.*.bam", recursive=T, include.dirs=T)
bam = bamfiles[grep("bam$", bamfiles)]
bai = bamfiles[grep("bai", bamfiles)]
sampleid = gsub("/.*", "", bam)
df = cbind.data.frame(sampleid, bam, bai)
df = df[df$sampleid %in% samples,]

# Generate directories for each sample
lapply(df[,1], function(x) dir.create(paste0(wd, x)))

# Vcf files
# Mutect
files <- list.files(path=mutect.folder, pattern="PASS.vcf", recursive=T, include.dirs=T, full.names = T)
sampleid = basename(dirname(files))
df1 = cbind.data.frame(sampleid, files)
df1 = df1[df1$sampleid %in% samples,]
# Strelka
files <- list.files(path=strelka.folder, pattern="passed.somatic.indels.vcf$", recursive=T, include.dirs=T, full.names = T)
sampleid = basename(dirname(files))
df2 = cbind.data.frame(sampleid, files)
df2 = df2[df2$sampleid %in% samples,]
# Combine vcfs per sample
dfsplit <- rbind(df1,df2)
dfsplit0 <- split(dfsplit, f = as.character(dfsplit[,1]))
summary(names(dfsplit0) == df[,1])


# Yaml file
ethnicity <- "Asian"
epitope <- "8,9,10,11"

for(i in 1:nrow(df))
{
  id.line <- paste0('  ',df[i,1],':')
  bam.line <- paste0(bam.folder, df[i,2])
  bai.line <- paste0(bam.folder,df[i,3])
  
  writeLines(paste0('ethnicity: ',ethnicity),paste0(wd,df[i,1],'/params.yaml'))
  write(paste0('epitope_length: ',epitope),paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
  write('samples:',paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
  write(id.line,paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
  write(paste0('    normal_bam: ', bam.line),paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
  write(paste0('    normal_bai: ', bai.line),paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
  write('    vcfs:',paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
  write(paste0('    - ',dfsplit0[[i]][,2]),paste0(wd,df[i,1],'/params.yaml'),append = TRUE)
}

