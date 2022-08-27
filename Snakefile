rule all:
    input:
        expand("out/{sample}/polysolver/winners.hla.txt",sample=config['samples']),
        expand("out/{sample}/polysolver/winners.hla.line.txt",sample=config['samples']),
        expand("out/{sample}/pvacseq/all.snp.flt.tsv",sample=config['samples']),
        expand("out/{sample}/pvacseq/all.indel.flt.tsv",sample=config['samples'])
        

rule link_bam:
    input:
        bam = lambda wildcards: config['samples'][wildcards.sample]['normal_bam'],
        bai = lambda wildcards: config['samples'][wildcards.sample]['normal_bai']
    output:
        bam = "out/{sample}/bam/normal.bam",
        bai = "out/{sample}/bam/normal.bam.bai"
    threads:
        2
    message:
        "Making soft-links to bam files"
    shell:
        "ln -s {input.bam} {output.bam};"
        "ln -s {input.bai} {output.bai}"


rule polysolver:
    input:
        bam = "out/{sample}/bam/normal.bam"
    params:
        ethnicity = config['ethnicity']
    output:
        winners = "out/{sample}/polysolver/winners.hla.txt",
        winners_line = "out/{sample}/polysolver/winners.hla.line.txt"
    log:
        "out/{sample}/polysolver/polysolver.log"
    threads:
        8
    message:
        "Running HLA typing"
    shell:
        "source /mnt/projects/phuazjc/workspace/software/polysolver/scripts/config.bash;"
        "/mnt/projects/phuazjc/workspace/software/polysolver/scripts/shell_call_hla_type {input.bam} {params.ethnicity} 1 hg19 STDFQ 0 out/{wildcards.sample}/polysolver;"
        """rm out/{wildcards.sample}/polysolver/*temp*;"""
        """rm out/{wildcards.sample}/polysolver/*lik*;"""
        """rm out/{wildcards.sample}/polysolver/*R0k6*;"""
        """awk '{{print $2","$3}}' {output.winners} | awk '{{ print toupper($0) }}' | awk -F',' -v OFS=',' '{{split($1,a,"_"); print a[1]"-"a[2]"*"a[3]":"a[4]}}' > out/{wildcards.sample}/polysolver/allelle1;"""
        """awk '{{print $2","$3}}' {output.winners} | awk '{{ print toupper($0) }}' | awk -F',' -v OFS=',' '{{split($2,a,"_"); print a[1]"-"a[2]"*"a[3]":"a[4]}}' > out/{wildcards.sample}/polysolver/allelle2;"""
        """paste out/{wildcards.sample}/polysolver/allelle1 out/{wildcards.sample}/polysolver/allelle2 > out/{wildcards.sample}/polysolver/alleles;"""
        """awk -vORS=' ' 1 out/{wildcards.sample}/polysolver/alleles | awk '{{for(i=1;i<=NF;i++) a[$i]++}} END{{for(i in a) printf i" ";print ""}}' | tr ' ' ',' | sed 's/,$//' > {output.winners_line};"""
        """rm out/{wildcards.sample}/polysolver/allele*"""

#############################
### MODIFIED FROM HECHUAN ###
#############################

rule concat_vcf:
    input:
        lambda wildcards: config['samples'][wildcards.sample]['vcfs']
    output:
        "out/{sample}/pvacseq/{sample}.concat.vcf"
    threads:
        2
    message:
        "Concatenating vcfs"
    shell:
        '''cat {input} |awk '!/^#/ && $1~/^(chr)?[1-9XYM]/ && $7=="PASS" ' '''
        ''' |cut -f1-7|sort -k1,1 -k2,2n |uniq '''
        ''' |awk 'BEGIN{{print "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER"}}1' > {output}'''


rule vep:
    input:
        "out/{sample}/pvacseq/{sample}.concat.vcf"
    output:
        "out/{sample}/pvacseq/{sample}.vep.vcf"
    log:
        "out/{sample}/pvacseq/{sample}.vep.log"
    threads:
        8
    message:
        "Annotating variants"
    shell:
        "/home/phuazjc/bin/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl --input_file {input} --format vcf --output_file {output} --vcf"
        " --symbol --terms SO --plugin Downstream --plugin Wildtype --offline --pick --force_overwrite > {log}"


rule pvacseq:
    params:
        epitope_length=lambda wildcards: config['epitope_length']
    input:
        "out/{sample}/pvacseq/{sample}.vep.vcf"
    output:
        "out/{sample}/pvacseq/MHC_Class_I/{sample}.final.tsv"
    log:
        "out/{sample}/pvacseq/{sample}.pvac.log"
    threads:
        8
    message:
        "Running neoantigen prediction"
    shell:
        """allele=$(cat out/{wildcards.sample}/polysolver/winners.hla.line.txt);"""
        "pvacseq run -e {params.epitope_length} --iedb-install-directory /mnt/projects/yangh/tcga/bin/software_download/pvac/"
        " {input} {wildcards.sample} $allele NetMHCpan out/{wildcards.sample}/pvacseq/ > {log}"


rule binding_filter:
    input:
        "out/{sample}/pvacseq/MHC_Class_I/{sample}.final.tsv"
    output:
        "out/{sample}/pvacseq/MHC_Class_I/{sample}.final.flt.tsv"
    threads:
        8
    message:
        "Keeping neoantigens with <500nM binding affinity"
    shell:
        """perl -lane 'next if $.==1 || $F[18]>=500; print join("\\t","{wildcards.sample}",@F[0,1,3,4,7..16,18,19])' {input}"""
        """ |sort -k2,2 -k3,3n > {output}"""


rule all_snp:
    input:
        expand('out/{sample}/pvacseq/MHC_Class_I/{sample}.final.flt.tsv',sample=config['samples'].keys())
    output:
        'out/{sample}/pvacseq/all.snp.flt.tsv'
    threads:
        8
    message:
        "Outputting neoantigens generated from SNPs"
    shell:
        """awk 'BEGIN{{print "Sample\\tChromosome\\tPosition\\tReference\\tVariant\\tVariantType\\tMutation\\tProteinPosition\\tGeneName\\tHLA_Allele\\tPeptideLength\\tSub-peptide_Position\\tMutationPosition\\tMT_Epitope_Seq\\tWT_Epitope_Seq\\tNetMHCpan_MT_Score\\tCorresponding_WT_Score"}} length($4)==1 && length($5)==1' {input} > {output}"""


rule all_indel:
    input:
        expand('out/{sample}/pvacseq/MHC_Class_I/{sample}.final.flt.tsv',sample=config['samples'].keys())
    output:
        'out/{sample}/pvacseq/all.indel.flt.tsv'
    threads:
        8
    message:
        "Outputting neoantigens generated from indels"
    shell:
        """awk 'BEGIN{{print "Sample\\tChromosome\\tPosition\\tReference\\tVariant\\tVariantType\\tMutation\\tProteinPosition\\tGeneName\\tHLA_Allele\\tPeptideLength\\tSub-peptide_Position\\tMutationPosition\\tMT_Epitope_Seq\\tWT_Epitope_Seq\\tNetMHCpan_MT_Score\\tCorresponding_WT_Score"}} length($4)!=1 || length($5)!=1' {input} > {output}"""

