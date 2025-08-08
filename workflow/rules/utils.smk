rule samtools_index:
    input:
        "{file_path}.bam",
    output:
        "{file_path}.bai",
    log:
        "logs/samtools_index/{file_path}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v7.2.0/bio/samtools/index"
