from Bio import SeqIO

#samples list
SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

#the final goal is to generate the differential expression results and the full pipeline report
rule all:
    input:
        "GCF_000845245.1_results.txt",
        "PipelineReport.txt"

#Rule 1: getting the HCMV genome CDS from NCBI
rule getting_cds:
    output:
        zip = "GCF_000845245.1.zip",
        raw_fasta = "raw_cds.fna"
    shell:
        """
        datasets download genome accession GCF_000845245.1 --include cds --filename {output.zip}
        unzip -p {output.zip} ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna > {output.raw_fasta}
        """

# Rule 2: clean the FASTA and count how many CDS are in the genome
rule format_file:
    input: "raw_cds.fna"
    output:
        clean = "GCF_000845245.1_cds_clean.fasta",
        report = "initial_cds_count.txt"
    run:
        count = 0
        with open(output.clean, "w") as f:
            for rec in SeqIO.parse(input[0], "fasta"):
                if "[protein_id=" in rec.description:
                    pid = rec.description.split("[protein_id=")[1].split("]")[0]
                    f.write(f">{pid}\n{rec.seq}\n")
                    count += 1
        with open(output.report, "w") as f:
            f.write(f"The HCMV genome (GCF_000845245.1) has {count} CDS.\n")

# Rule 3: create kallisto index from the cleaned CDS and quantification of TPM
rule kallisto:
    input: "GCF_000845245.1_cds_clean.fasta"
    output: "GCF_000845245.1_index.idx"
    shell: "kallisto index -i {output} {input}"

rule kallisto_quant:
    input: 
        index = "GCF_000845245.1_index.idx",
        r1 = "data/{sample}_1.fastq",
        r2 = "data/{sample}_2.fastq"
    output: directory("results/{sample}")
    shell: "kallisto quant -i {input.index} -o {output} -b 30 {input.r1} {input.r2}"

# Rule 4: running sleuth in R to find differentially expressed genes (FDR < 0.05)
rule sleuth:
    input:
        meta = "metadata.txt",
        kallisto = expand("results/{sample}", sample=SAMPLES)
    output: "GCF_000845245.1_results.txt"
    shell: "Rscript run_sleuth.r {output}"

#Rule 5: build a Bowtie2 index to filter reads and map reads (save only the mapped pairs)
rule bowtie_index:
    input: "raw_cds.fna"
    output: "hcmv_idx.1.bt2"
    shell: "bowtie2-build {input} hcmv_idx"

rule bowtie_map:
    input:
        r1 = "data/{sample}_1.fastq",
        r2 = "data/{sample}_2.fastq",
        index = "hcmv_idx.1.bt2"
    output:
        sam = "mapped_{sample}.sam",
        f1 = "mapped_{sample}_1.fq",
        f2 = "mapped_{sample}_2.fq"
    shell: "bowtie2 -x hcmv_idx -1 {input.r1} -2 {input.r2} -S {output.sam} --al-conc mapped_{wildcards.sample}_%.fq"

#Rule 6: count read pairs before and after filtering (mapping statistics)
rule map_stats:
    input:
        before = "data/{sample}_1.fastq",
        after = "mapped_{sample}_1.fq"
    output: 
        stats = "{sample}_mapping_stats.txt"
    run:
        bef = sum(1 for l in open(input.before)) / 4
        aft = sum(1 for l in open(input.after)) / 4
        with open(output.stats, "w") as f:
            f.write(f"Sample {wildcards.sample} had {int(bef)} read pairs before and {int(aft)} after Bowtie2 filtering.\n")

#Rule 7: assemble filtered reads into contigs using SPAdes (--rna for transcriptomes)
rule spades:
    input:
        r1 = "mapped_{sample}_1.fq",
        r2 = "mapped_{sample}_2.fq"
    output: "assembly_{sample}/transcripts.fasta"
    shell: "spades.py -k 127 --rna -1 {input.r1} -2 {input.r2} -o assembly_{wildcards.sample}"

#Rule 8: find the longest contig from the assembly using a loop
rule longest_contig:
    input: "assembly_{sample}/transcripts.fasta"
    output: "longest_contig_{sample}.fasta"
    run:
        #load all sequences into a list
        records = list(SeqIO.parse(input[0], "fasta"))
        
        if records:
            #start by picking the first one
            longest_rec = records[0]
            
            #look at every sequence one by one
            for rec in records:
                #if this one is longer than current
                if len(rec.seq) > len(longest_rec.seq):
                    #then it replaces the old one
                    longest_rec = rec
            
            #save only the last one
            SeqIO.write(longest_rec, output[0], "fasta")

# Rule 9: creating a database for blast
rule blast_db:
    output: "Betaherpesvirinae_db.nhr"
    shell:
        """
        datasets download virus genome taxon 10357 --filename virus.zip
        unzip -p virus.zip "*.fna" > Betaherpesvirinae.fasta
        makeblastdb -in Betaherpesvirinae.fasta -dbtype nucl -out Betaherpesvirinae_db
        """

#Rule 10: blast the longest contig against the database
rule blast:
    input:
        query = "longest_contig_{sample}.fasta",
        db = "Betaherpesvirinae_db.nhr" 
    output: "{sample}_blast.txt"
    shell:
        """
        #sample header
        echo "{wildcards.sample}:" > {output}

        #column headers
        echo -e "sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle" >> {output}

        #run blast command and append the results
        blastn -query {input.query} -db Betaherpesvirinae_db \
        -max_target_seqs 5 -max_hsps 1 \
        -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' \
        >> {output}
        """

#Rule 11: writing final report with the results 
rule final_report:      
    input:
        cds = "initial_cds_count.txt",
        sleuth = "GCF_000845245.1_results.txt",
        stats = expand("{sample}_mapping_stats.txt", sample=SAMPLES),
        blast = expand("{sample}_blast.txt", sample=SAMPLES)
    output: "PipelineReport.txt"
    run:
        with open(output[0], "w") as out:
            #CDS count
            out.write(open(input.cds).read() + "\n")
            
            #sleuth
            out.write("target_id\ttest_stat\tpval\tqval\n")
            s_lines = open(input.sleuth).readlines()
            if len(s_lines) > 1: 
                out.write("".join(s_lines[1:]))
            
            #mapping stats
            for s in input.stats: 
                out.write(open(s).read())
            
            #blast results
            for b in input.blast: 
                out.write(open(b).read() + "\n")