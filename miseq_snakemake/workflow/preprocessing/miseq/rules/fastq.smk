localrules: bcl_to_fastq, aggregater,remove_ill_id 

def ill_id_removal(wildcards):
    """creates the input for the ampliseq_quant function while also removing the 
    ill_id from the file names"""
    miseq_dir = samp_dict[wildcards.sample]
    samp_sheet = samp_sheet_dict[miseq_dir]
    if config['bad_run']:
        return 'data/aligned/Undetermined_S0.bam'

    else:
        sub_samps = read_sample_sheet(samp_sheet)['Sample_ID'].tolist()
        ill_id = sub_samps.index(wildcards.sample) + 1
        return expand('data/raw/fastq/{{sample}}_S{ill_id}_L001_R{{read}}_001.fastq.gz',ill_id=ill_id)[0]

rule bcl_to_fastq:
    input:
        sample_sheet = lambda wildcards: samp_sheet_dict[wildcards.miseq_dir],
        miseq_loc = os.path.join(config["miseq_output"],'{miseq_dir}')
    output:
        out_dir=directory('data/raw/{miseq_dir}/')
    conda:
        'envs/bcl2fastq2.yaml'
    params:
        get_flags(config['bcl2fastq_flags'])
    shell:
        """ 
        # if [ -f "{input.miseq_loc}/GenerateFASTQRunStatistics.xml" ] ; then
        #     echo "Fastq detected"
        #     mkdir {output.out_dir}
        #     ln -s {input.miseq_loc}/Alignment*/*/Fastq/*.fastq.gz {output.out_dir}/fastq
        # else 
            echo 'bcl2fastq -R {input.miseq_loc} --input-dir {input.miseq_loc}/Data/Intensities/BaseCalls/ --output-dir data/raw/{wildcards.miseq_dir}/ --sample-sheet {input.sample_sheet} {params}';
            bcl2fastq -R {input.miseq_loc} --input-dir {input.miseq_loc}/Data/Intensities/BaseCalls/ --output-dir data/raw/{wildcards.miseq_dir}/ --sample-sheet {input.sample_sheet} {params}
        # fi
        """

rule aggregater:
    input:
        miseq_dirs=expand('data/raw/{miseq_dir}/',miseq_dir=miseq_dirs)
    output:
        read="data/raw/fastq/{sample}_S{ill_id}_L001_R{read}_001.fastq.gz"
    shell:
        # lazy implementation, fix when have time
        """
        sample={wildcards.sample}_S{wildcards.ill_id}_L001_R{wildcards.read}_001.fastq.gz
        fastq_path=$(find -name $sample)
        ln -rs $fastq_path {output.read}
        """

rule remove_ill_id:
    input:
        bam_file = ill_id_removal
    output:
        no_ill_id='data/raw/fastq/{sample}_R{read}.fastq.gz'
    shell:
        "mv {input} {output}" 