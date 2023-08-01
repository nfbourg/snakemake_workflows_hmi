localrules: fastqc, multiqc

def multiqc_input(bad_run):
    """creates the input for the ampliseq_quant function while also removing the 
    ill_id from the file names"""
    ill_ids = []
    for sample in samps:
        miseq_dir = samp_dict[sample]
        samp_sheet = samp_sheet_dict[miseq_dir]
        if bad_run:
            sub_samps = ['Undetermined']
            ill_ids.append(sub_samps.index(sample))

        else:
            sub_samps = read_sample_sheet(samp_sheet)['Sample_ID'].tolist() 
            ill_ids.append(sub_samps.index(sample) + 1)

    return expand(expand("qc/fastqc/{sample}_R{{read}}_fastqc.zip",zip,
                          sample=samps,ill_id=ill_ids), read=reads)

rule fastqc:
    input:
        "data/raw/fastq/{sample}_R{read}.fastq.gz"
    output:
        html="qc/fastqc/{sample}_R{read}.html",
        zip="qc/fastqc/{sample}_R{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}_R{read}.log"
    threads: 1
    wrapper:
        "v1.5.0/bio/fastqc"

rule multiqc:
    input:
        multiqc_input(config['bad_run'])
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"